#include "FunctionTable.h"
#include "RbFileManager.h"
#include "RevClient.h"
#include "RlFunction.h"
#include "RlUserInterface.h"
#include "Parser.h"
#include "Workspace.h"
#include "RbSettings.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Environment.h"
#include "RevPtr.h"
#include "TypeSpec.h"
#include "boost/algorithm/string/trim.hpp"
#include <unistd.h>

#ifdef RB_MPI
#include <mpi.h>
#endif

#include "replxx.hxx"

#include <ranges>
#include <filesystem>
#include <boost/algorithm/string/predicate.hpp>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "RbException.h"
//#define ctrl(C) ((C) - '@')

const char* default_prompt = (char*) "> ";
const char* incomplete_prompt = (char*) "+ ";
const char* prompt = default_prompt;

using namespace RevLanguage;

namespace fs = std::filesystem;

std::vector<std::string> getFileList(const RevBayesCore::path& dir)
{
    std::vector<RevBayesCore::path> filenames;

    RevBayesCore::setStringWithNamesOfFilesInDirectory( RevBayesCore::current_path() / dir, filenames, false );

    std::vector<std::string> v;

    for(auto& filename: filenames)
        v.push_back( filename.string() );
    
    return v;
}



std::vector<std::string> getDefaultCompletions( void )
{
    using namespace RevLanguage; 

    std::set<std::string> all_names;
    
    const FunctionTable& ft = RevLanguage::Workspace::userWorkspace().getFunctionTable();
    for (auto& [func_name, func]: ft)
    {
        all_names.insert(func_name);
    }
    
    // This seems to be a duplicate??
    std::vector<std::string> function_table_names;
    ft.getFunctionNames(function_table_names);
    for (auto& func_name: function_table_names)
    {
        all_names.insert(func_name);
    }
    
    for (auto& [var_name, var]: Workspace::userWorkspace().getVariableTable())
    {
        all_names.insert(var_name);
    }
    
    for (auto& [var_name, var]: Workspace::globalWorkspace().getVariableTable())
    {
        all_names.insert(var_name);
    }

    for (auto& [type_name, type]: Workspace::userWorkspace().getTypeTable())
    {
        all_names.insert(type_name);
    }
    
    for (auto& [type_name, type]: Workspace::globalWorkspace().getTypeTable())
    {
        all_names.insert(type_name);
    }
    
    return all_names | std::ranges::to<std::vector>();
}

replxx::Replxx::completions_t getCompletions(const std::string& buf, int& contextLen)
{
    std::string cmd = buf;
    std::vector<std::string> completions;

    // parse command
    RevLanguage::ParserInfo pi =
        RevLanguage::Parser::getParser().checkCommand(
            cmd,
            RevLanguage::Workspace::userWorkspacePtr()
            );

    if (pi.inComment) {
        return replxx::Replxx::completions_t{};
    }

    size_t commandPos = 0;

    if (pi.inQuote) {
        // ---------- in quote ------------
        commandPos = cmd.rfind("\"") + 1;
        completions = getFileList(cmd.substr(commandPos));
    }
    else {
        std::vector<std::string> expressionSeparator = {
            " ", "%", "~", "=", "&", "|", "+", "-", "*", "/",
            "^", "!", "=", ",", "<", ">", ")", "[", "]", "{", "}"
        };

        for (auto const& s : expressionSeparator) {
            size_t find_idx = cmd.rfind(s);
            if (find_idx < cmd.size()) {
                commandPos = std::max(commandPos, find_idx);
            }
        }

        if (pi.function_name != "") {
            if (pi.argument_label != "") {
                commandPos = cmd.rfind("=") + 1;
                completions = getDefaultCompletions();
            }
            else {
                commandPos = std::max(cmd.rfind("("), cmd.rfind(",")) + 1;

                auto v = Workspace::globalWorkspace()
                    .getFunctionTable()
                    .findFunctions(pi.function_name);

                for (auto* fn : v) {
                    const auto& argRules = fn->getArgumentRules();
                    for (size_t i = 0; i < argRules.size(); i++) {
                        completions.push_back(argRules[i].getArgumentLabel());
                    }
                }
            }
        }
        else {
            if (commandPos > 0) {
                commandPos++;
            }
            completions = getDefaultCompletions();
        }
    }

    // skip leading spaces
    while (commandPos < cmd.size() && cmd[commandPos] == ' ') {
        commandPos++;
    }

    // 🔑 THIS replaces your prefix logic
    contextLen = static_cast<int>(cmd.size() - commandPos);

    std::string compMatch = cmd.substr(commandPos);

    // No completions if we don't have at least two characters.
    if (compMatch.size() < 2) return {};

    replxx::Replxx::completions_t result;

    for (auto const& m : completions) {
        if (boost::starts_with(m, compMatch)) {
            result.emplace_back(m);
        }
    }

    return result;
}

replxx::Replxx::hints_t getHints(std::string const& buf,
                                 int& contextLen,
                                 replxx::Replxx::Color& color)
{
    std::string cmd = buf;

    RevLanguage::ParserInfo pi =
        RevLanguage::Parser::getParser().checkCommand(
            cmd,
            RevLanguage::Workspace::userWorkspacePtr()
            );

    replxx::Replxx::hints_t hints;

    if (!pi.inQuote &&
        Workspace::globalWorkspace().existsFunction(pi.function_name))
    {
        auto functions =
            Workspace::globalWorkspace()
            .getFunctionTable()
            .findFunctions(pi.function_name);

        for (auto* f : functions)
        {
            std::string hint;

            hint += f->getReturnType().getType();
            hint += " ";
            hint += pi.function_name;
            hint += "(";

            const auto& argRules = f->getArgumentRules();
            for (size_t i = 0; i < argRules.size(); i++)
            {
                hint += argRules[i].getArgumentLabel();
                if (i + 1 < argRules.size())
                    hint += ", ";
            }

            hint += ")";

            hints.emplace_back(std::move(hint));
        }

        contextLen = 0;
        color = replxx::Replxx::Color::GRAY;
    }

    return hints;
}

/**
 * tab completion callback
 * 
 * Update list of available completions.
 * 
 * @param buf
 * @param lc

void completeOnTab(const char *buf, linenoiseCompletions *lc)
{
    bool debug = false;
    std::string cmd = buf;
    std::vector<std::string> completions;

    // parse command
    RevLanguage::ParserInfo pi = RevLanguage::Parser::getParser().checkCommand(cmd, RevLanguage::Workspace::userWorkspacePtr());

    if (pi.inComment)
    {
        if (debug) { std::cout << "linenoise-debug: pi.inComment==TRUE\n"; }
        // no completions available in comments
        return;
    }

    // set completions and position on command line where to start matching completions
    size_t commandPos = 0;
    if (pi.inQuote)
    {
        if (debug) { std::cout << "linenoise-debug: pi.inQuote==TRUE\n"; }
        // ---------- in quote ------------
        // search for files with portion after the opening quote                
        commandPos = cmd.rfind("\"") + 1;
        completions = getFileList(cmd.substr(commandPos, cmd.size()));
    }
    else
    {
        if (debug) { std::cout << "linenoise-debug: pi.inComment==FALSE\n"; }
        std::vector<std::string> expressionSeparator;
        
        expressionSeparator.push_back(" ");
        expressionSeparator.push_back("%");
        expressionSeparator.push_back("~");
        expressionSeparator.push_back("=");
        expressionSeparator.push_back("&");
        expressionSeparator.push_back("|");
        expressionSeparator.push_back("+");
        expressionSeparator.push_back("-");
        expressionSeparator.push_back("*");
        expressionSeparator.push_back("/");
        expressionSeparator.push_back("^");
        expressionSeparator.push_back("!");
        expressionSeparator.push_back("=");
        expressionSeparator.push_back(",");
        expressionSeparator.push_back("<");
        expressionSeparator.push_back(">");
        expressionSeparator.push_back(")");
        expressionSeparator.push_back("[");
        expressionSeparator.push_back("]");
        expressionSeparator.push_back("{");
        expressionSeparator.push_back("}");

        // find position of right most expression separator in cmd

        for(auto& s: expressionSeparator)
        {
            if (debug) { std::cout << "linenoise-debug: rfind(\"" << s << "\",\"" << cmd << "\")=" << cmd.rfind(s) << "\n"; }
            size_t find_idx = cmd.rfind(s);
            if (find_idx < cmd.size())
                commandPos = std::max(commandPos, cmd.rfind(s));
        }
        if (debug) { std::cout << "linenoise-debug: cmd.size()=" << cmd.size() << "\n"; }
        if (debug) { std::cout << "linenoise-debug: commandPos=" << commandPos << "\n"; }

        // special hack: for some reason, base_variable is only set by the parser when there is no trailing characters after the dot
        // find position of right most dot
        // Sebastian: Currently unused
//        size_t dotPosition = cmd.rfind(".");

        if (pi.function_name != "")
        {
            if (debug) { std::cout << "linenoise-debug: pi.function_name!=\"\"\n"; }
            // ---------- function defined ------------
            if (pi.argument_label != "") // assigning an argument label
            {
                if (debug) { std::cout << "linenoise-debug: pi.argument_label!=\"\"\n"; }
                commandPos = cmd.rfind("=") + 1;
                // not sure exactly what should be here... setting completions to everything
                completions = getDefaultCompletions();

            }
            else // break on either '(' or ','
            {
                if (debug) { std::cout << "linenoise-debug: pi.argument_label==\"\"\n"; }
                commandPos = std::max(cmd.rfind("("), cmd.rfind(",")) + 1;
                
                std::vector<Function *> v = Workspace::globalWorkspace().getFunctionTable().findFunctions(pi.function_name);
                
                for (std::vector<Function *>::iterator it = v.begin(); it != v.end(); it++)
                {
                    const RevLanguage::ArgumentRules& argRules = (*it)->getArgumentRules();
                    for (size_t i = 0; i < argRules.size(); i++)
                    {
                        completions.push_back(argRules[i].getArgumentLabel());
                    }
                }
            }
        }
        else
        {
            if (debug) { std::cout << "linenoise-debug: pi.function_name==\"\"\n"; }
            // ---------- default -----------            
            if (commandPos > 0)
            {
                commandPos++;
            }
            completions = getDefaultCompletions();

        }
    }

    // discard any extra space in beginning of the string that is used to match against completions
    while (buf[commandPos] == ' ')
    {
        commandPos++;
    }

    // match partial command and pass filtered completions to linenoise
    std::string previousCommands;
    for (int i = 0; i < commandPos; i++)
    {
        previousCommands += buf[i];
    }

    // the string the matching is made against
    std::string compMatch(buf + commandPos);

    // populate linenoise completions
    
    for(auto& m: completions)
    {
        if (boost::starts_with(m, compMatch))
        {
            std::string c = previousCommands + m;
            linenoiseAddCompletion(lc, c.c_str());
        }
    }
}

int printFunctionParameters(const char *buf, size_t len, char c)
{
    std::string cmd = buf;
    RevLanguage::ParserInfo pi = RevLanguage::Parser::getParser().checkCommand(cmd, RevLanguage::Workspace::userWorkspacePtr());
    if ( Workspace::globalWorkspace().existsFunction(pi.function_name) )
    {

        std::vector<Function *> functions = Workspace::globalWorkspace().getFunctionTable().findFunctions(pi.function_name);
        
        for (std::vector<Function *>::iterator it = functions.begin(); it != functions.end(); ++it)
        {
            Function *f = *it;
            std::cout << "\n\r" + f->getReturnType().getType() + " " + pi.function_name + " (";

            const RevLanguage::ArgumentRules& argRules = (*it)->getArgumentRules();
            for (size_t i = 0; i < argRules.size(); i++)
            {
                std::cout << argRules[i].getArgumentLabel();
                if (i < argRules.size() - 1) {
                    std::cout << ", ";
                }
            }

        }
        
        
        std::cout << ")\n\r";
        linenoiceSetCursorPos(0);
        std::cout << prompt << buf;
        std::cout.flush();
    }
    return 0;
}

*/

namespace RevClient
{

int interpret(const std::string& command)
{
    size_t bsz = command.size();
#ifdef RB_MPI
    MPI_Bcast(&bsz, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    char * buffer = new char[bsz+1];
    buffer[bsz] = 0;
    for (int i = 0; i < bsz; i++)
        buffer[i] = command[i];
#ifdef RB_MPI
    MPI_Bcast(buffer, (int)bsz, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

    std::string tmp = std::string( buffer );

    return RevLanguage::Parser::getParser().processCommand(tmp, RevLanguage::Workspace::userWorkspacePtr());
}

void shutdown()
{
    Workspace::userWorkspace().clear();
    Workspace::globalWorkspace().clear();

#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
}


void execute_file(const fs::path& filename, bool echo, bool continue_on_error)
{
    auto& settings = RbSettings::userSettings();
    std::stringstream inFile = RevBayesCore::readFileAsStringStream(filename);

    // Command-processing loop
    std::string commandLine;
    int lineNumber = 0;
    int result = 0;     // result from processing of last command
    while ( inFile.good() )
    {
        // Read a line
        std::string line;
        RevBayesCore::safeGetline(inFile, line);
        lineNumber++;
        
        if ( echo )
        {
            if ( result == 1 )
            {
                std::cout << "+ " << line << std:: endl;
            }
            else
            {
                std::cout << "> " << line << std::endl;
            }
        }
        
        // If previous result was 1 (append to command), we do this
        if ( result == 1 )
        {
            commandLine += line;
        }
        else
        {
            commandLine = line;
        }

        // Process the line and record result
        result = Parser::getParser().processCommand( commandLine, Workspace::userWorkspacePtr() );
        if ( result == 2 )
        {
            if (not continue_on_error)
                throw RbException() << "Problem processing line " << lineNumber << " in file " << filename;
            else
            {
                std::ostringstream err;
                err<<"Error:\tProblem processing line " << lineNumber << " in file " << filename;
                RBOUT(err.str());
            }
        }
    }
}

/**
 * Main application loop.
 * 
 */
void startInterpreter()
{
    auto& settings = RbSettings::userSettings();

    // If we aren't using MPI, this will be zero.
    // If we are using MPI, it will be zero for the first process.
    int pid = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif

    replxx::Replxx rx;

    if ( pid == 0 )
    {
        rx.set_word_break_characters(" %~=&|+-*/^!,<>\t\n\"\\'`@$;{(");

        /* Load history from file. The history file is just a plain text file
         * where entries are separated by newlines. */
        rx.history_load("history.txt");

        /* Set tab completion callback */
        rx.set_completion_callback( getCompletions );

        /* Callback for printing function signatures on '('.
           It doesn't seem to do anything at the moment: pi.function_name is never set. */
        rx.set_hint_callback( getHints );

//        // Replxx uses CTRL-P/N to cycle through alternative completions.
//        // These settings make Replxx match linenoise.        
//
//        rx.bind_key_internal( replxx::Replxx::KEY::control('P'),
//                              "history_previous" );
//
//        rx.bind_key_internal( replxx::Replxx::KEY::control('N'),
//                              "history_next" );
    }

    int result = 0;
    std::string commandLine;
    std::string cmd;

    while (true)
    {
        
        // set prompt
        if (result == 0 || result == 2)
        {
            prompt = default_prompt;
        }
        else //if (result == 1)
        {
            prompt = incomplete_prompt;
        }

        // process command line
        if ( pid == 0 )
        {
            const char* line = rx.input("> ");

            if (!line) 
            {
                // [JS, 2015-11-03]
                // Null input, e.g. if the user entered CTRL-D or CTRL-C.
                // If not handled here, segmentation fault results. Not a dealbreaker, but annoying.
                shutdown();

                exit(0);
            }
            else if (*line)
            {
                rx.history_add(line);                /* Add to the history. */
                rx.history_save("history.txt");      /* Save the history on disk. */
           
                cmd = line;
                boost::trim(cmd);

                if (cmd == "clr" || cmd == "clear")
                {
                    rx.clear_screen();
                }
                else
                {
                    // interpret Rev statement
                    if (result == 0 || result == 2)
                    {
                        commandLine = cmd;
                    }
                    else //if (result == 1)
                    {
                        commandLine += "\n " + cmd;
                    }
                }
            }
        }
        
        result = interpret(commandLine);

        if (result == 2)
        {
            commandLine.clear();
        }

/* The typed string is returned as a malloc() allocated string by
         * linenoise, so the user needs to free() it. */
    }
    
    
}

void startJupyterInterpreter()
{
    // If we aren't using MPI, this will be zero.
    // If we are using MPI, it will be zero for the first process.
    int pid = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
    
    /* Declare things we need */
    int result = 0;
    std::string commandLine = "";

    for (;;)
    {
        /* Print prompt based on state after previous iteration */
        if ( pid == 0 )
        {
            if (result == 0 || result == 2)
            {
                std::cout << "> ";
            }
            else
            {
                std::cout << "+ ";
            }

            /* Get the line */
            std::string line = "";
            if (not std::getline(std::cin, line)) return;

            if (result == 0 || result == 2)
            {
                commandLine = line;
            }
            else if (result == 1)
            {
                commandLine += ";" + line;
            }
        }

        result = RevClient::interpret(commandLine);
    }
}

}
