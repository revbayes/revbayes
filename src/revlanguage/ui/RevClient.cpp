#include "FunctionTable.h"
#include "RbFileManager.h"
#include "RevClient.h"
#include "RlFunction.h"
#include "Parser.h"
#include "Workspace.h"
#include "RbSettings.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Environment.h"
#include "RevPtr.h"
#include "TypeSpec.h"
#include "boost/algorithm/string/trim.hpp"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifdef RB_MPI
#include <mpi.h>
#endif

extern "C" {
#include "linenoise.h"
}

#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
//#define ctrl(C) ((C) - '@')

const char* default_prompt = (char*) "> ";
const char* incomplete_prompt = (char*) "+ ";
const char* prompt = default_prompt;

using namespace RevLanguage;

namespace fs = boost::filesystem;

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
    std::set<std::string> c;
    
    const FunctionTable& ft = RevLanguage::Workspace::userWorkspace().getFunctionTable();
    for (std::multimap<std::string, Function*>::const_iterator it = ft.begin(); it != ft.end(); ++it)
    {
        c.insert(it->first);
    }
    
    std::vector<std::string> function_table_names;
    ft.getFunctionNames(function_table_names);
    for (size_t i = 0; i < function_table_names.size(); i++)
    {
//        std::cout << function_table_names[i] << "\n";
        c.insert(function_table_names[i]);
    }
    
    VariableTable v = RevLanguage::Workspace::userWorkspace().getVariableTable();
    
    for (VariableTable::iterator it = v.begin(); it != v.end(); ++it)
    {
        c.insert(it->first);
    }
    
    v = RevLanguage::Workspace::globalWorkspace().getVariableTable();
        
    for (VariableTable::iterator it = v.begin(); it != v.end(); ++it)
    {
        c.insert(it->first);
    }

    const TypeTable& t_user = RevLanguage::Workspace::userWorkspace().getTypeTable();
    
    for (TypeTable::const_iterator it = t_user.begin(); it != t_user.end(); ++it)
    {
        c.insert(it->first);
    }
    
    const TypeTable& t_global = RevLanguage::Workspace::globalWorkspace().getTypeTable();
    
    for (TypeTable::const_iterator it = t_global.begin(); it != t_global.end(); ++it)
    {
        c.insert(it->first);
    }
    
    std::vector<std::string> vec;
    for (std::set<std::string>::iterator it = c.begin(); it != c.end(); ++it)
    {
        vec.push_back( *it );
    }

    return vec;
}

/**
 * tab completion callback
 * 
 * Update list of available completions.
 * 
 * @param buf
 * @param lc
 */
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

        BOOST_FOREACH(std::string s, expressionSeparator)
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
    
    BOOST_FOREACH(std::string m, completions)
    {
        if (boost::starts_with(m, compMatch))
        {
            std::string c = previousCommands + m;
            linenoiseAddCompletion(lc, c.c_str());
        }
    }
}

/**
 * Print out function signatures
 * @param buf The buffer into which we print
 * @param len The length
 * @param c The character
 * @return  Returns the status
 */
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

/**
 * Main application loop.
 * 
 */
void startInterpreter( void )
{
    // If we aren't using MPI, this will be zero.
    // If we are using MPI, it will be zero for the first process.
    int pid = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
    
    /* Set tab completion callback */
    linenoiseSetCompletionCallback( completeOnTab );

    /* Load history from file. The history file is just a plain text file
     * where entries are separated by newlines. */
    if ( pid == 0 )
    {
        linenoiseHistoryLoad("history.txt"); /* Load the history at startup */
    }
    
    /* callback for printing function signatures on opening bracket*/

    // Currently disabled because
    // (i) it doesn't seem to do anything at the moment: pi.function_name is never set.
    // (ii) it makes parsing go crazy if the '(' is inside a string.

    // linenoiseSetCharacterCallback(printFunctionParameters, '(');

    int result = 0;
    std::string commandLine;
    std::string cmd;

    while (true)
    {
        
        char *line = NULL;
        
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
            line = linenoise(prompt);
            
            if (!line) 
            {
                // [JS, 2015-11-03]
                // Null input, e.g. if the user entered CTRL-D or CTRL-C.
                // If not handled here, segmentation fault results. Not a dealbreaker, but annoying.
                // There might be a more elegant way to do this, but right now, until I get
                // more familiar with the codebase or somebody else cares to improve it, we are
                // just going to replace the input with the intention, i.e. to quit, and let
                // the existing logic handle it.
                commandLine = "q();";
            }
            else
            {
                linenoiseHistoryAdd(line);              /* Add to the history. */
                linenoiseHistorySave("history.txt");    /* Save the history on disk. */
           
                cmd = line;
                boost::trim(cmd);

                if (cmd == "clr" || cmd == "clear")
                {
                    linenoiseClearScreen();
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

        /* The typed string is returned as a malloc() allocated string by
         * linenoise, so the user needs to free() it. */
        
        if ( pid == 0 )
        {
            free(line);
        }
        
    }
    
    
}

void startJupyterInterpreter( void )
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
