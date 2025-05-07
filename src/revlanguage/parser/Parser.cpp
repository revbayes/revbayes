#include "Parser.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <sstream>
#include <iostream>
#include <vector>

#include "RbException.h"
#include "RbHelpRenderer.h"
#include "RbHelpSystem.h"
#include "RbSettings.h"
#include "RbUtil.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "SyntaxElement.h"
#include "SyntaxFunctionCall.h"
#include "SyntaxVariable.h"
#include "RlUserInterface.h"
#include "Signals.h"
#include "Workspace.h"
#include "FunctionTable.h"
#include "RevObject.h"
#include "RlFunction.h"
#include "SyntaxFormal.h" // IWYU pragma: keep

#ifdef RB_MPI
#include <mpi.h>
#endif

#include "grammar.tab.h"

namespace RevBayesCore { class RbHelpEntry; }
namespace RevLanguage { class Environment; }

// Global flags indicating flex state
bool foundNewline;
bool foundErrorBeforeEnd;
bool foundEOF;


// Global command buffer RlString stream
std::stringstream rrcommand;


RevLanguage::Environment *executionEnvironment;

/** Constructor. Here we set the parser mode to executing. */
RevLanguage::Parser::Parser(void) {
    setParserMode(EXECUTING);
}

/**
 * This function breaks a buffer into Rev lines, taking into account the fact that ';' also
 * breaks lines in addition to newlines. The function also deals with different types of line
 * endings and translates all to Unix '\n' line endings.
 */
RevLanguage::ParserInfo RevLanguage::Parser::breakIntoLines(const std::string& cmd, std::vector<std::string>& lines) const {
    return breakIntoLines(cmd, lines, true);
}

RevLanguage::ParserInfo RevLanguage::Parser::breakIntoLines(const std::string& cmd, std::vector<std::string>& lines, bool validate) const {

    bool inComment, inQuote;

    // Initialize
    lines.clear(); // just in case
    inComment = inQuote = false;
    std::stringstream buf(cmd); // a stream is convenient for reading

    // Process command buffer
    while (buf.good()) {

        std::stringstream temp;
        bool escaped = false;

        while (buf.good()) {

            char c = char( buf.get());

            if (c == EOF && inQuote == true) {
                if (validate) {
                    throw RbException("End of line while in quote");
                }
                continue;
            } else if (escaped == false && c == '"') {
                /* switch quote on or off if not in comment */
                if (inQuote == true)
                    inQuote = false;
                else if (inComment == false)
                    inQuote = true;
            } else if (c == '#' && inQuote == false) {
                /* we are now in comment */
                inComment = true;
            } else if (c == ';' && inQuote == false && inComment == false) {
                /* break line here */
                break;
            } else if (c == EOF && inQuote == false) {
                /* break line here */
                break;
            } else if (c == '\n' && inQuote == false) {
                /* break line here */
                break;
            } else if (c == '\377' && inQuote == false) {
                /* break line here */
                //                c = '\n';
                continue;
            } else if (c == '\r' && inQuote == false) {
                /* break line here but first swallow any extra newline in DOS line ending */
                char d = char( buf.peek());
                if (d == '\n')
                    c = char( buf.get());
                break;
            }
            escaped = (c == '\\');

            temp.put(c);
        }


        if (temp.str().size() > 0)
        {
            lines.push_back(temp.str() + "\n");
        }
    }
    ParserInfo pi;
    pi.inComment = inComment;
    pi.inQuote = inQuote;
    pi.lines = lines;
    return pi;
}

/**
 * This function causes recursive execution of a syntax tree by calling the root to get its value.
 * As long as we return to the bison code, bison takes care of deleting the syntax tree. However,
 * if we encounter a quit() call, we delete the syntax tree ourselves and exit immediately.
 */
int RevLanguage::Parser::execute(SyntaxElement* root, Environment &env) const {

    // don't execute command if we are in checking mode
    if (RevLanguage::Parser::getParser().isChecking())
    {
        return 0;
    }

    // Declare a variable for the result
    RevPtr<RevVariable> result = NULL;

    //! Execute syntax tree
    try
    {
        result = root->evaluateContent(env);
    }
    catch (RbException& rbException)
    {

        std::ostringstream msg;

        // Catch a quit request
        if (rbException.getExceptionType() == RbException::QUIT)
        {
            delete( root);
            
            Workspace::userWorkspace().clear();
            Workspace::globalWorkspace().clear();
            
#ifdef RB_MPI
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
#endif
            
            exit(0);
        }

        // Catch a missing variable exception that might be interpreted as a request for
        // usage help on a function
        SyntaxVariable* rootPtr = dynamic_cast<SyntaxVariable*> ((SyntaxElement*) root);
        SyntaxVariable* theVariable = rootPtr;
        if (rbException.getExceptionType() == RbException::MISSING_VARIABLE && theVariable != NULL)
        {

            const std::string& fxnName = theVariable->getIdentifier();
            const std::vector<Function*>& functions = Workspace::userWorkspace().getFunctionTable().findFunctions(fxnName);
            if (functions.size() != 0)
            {
                for (std::vector<Function*>::const_iterator i = functions.begin(); i != functions.end(); i++)
                {
                    std::ostringstream s;
                    (*i)->printValue(s,true);
                    RBOUT(s.str());

                    // Uncommenting this as the function callSignature() does not produce the call signature despite its name
                    // -- Fredrik
                    // RBOUT( (*i)->callSignature() );
                }
                
                return 0;
            }
            
        }

        // All other exceptions
        rbException.print(msg);
        RBOUT(msg.str());

        // Return signal indicating problem
        return 2;
    }
    // Don't crash the intepreter for non-Rb exceptions such as file-not-found.
    catch (std::exception& e)
    {
        RBOUT(e.what());
        // Return signal indicating problem
        return 2;
    }

    // Print result if the root is not an assign expression
    if ( root->isAssignment() == false  &&  result != NULL  &&  result->getRevObject() != RevNullObject::getInstance() )
    {
        std::ostringstream msg;
        result->getRevObject().printValue(msg,true);
        RBOUT(msg.str());
    }

    // Warn if a return signal has been encountered
    // @todo Find out why the following lines do not work; they should
    //    if ( Signals::getSignals().isSet( Signals::RETURN ) )
    //        RBOUT( "WARNING: No function to return from" );
    Signals::getSignals().clearFlags();

    // Return success
    return 0;
}

/** Execute base variable expression to get the corresponding base variable */
void RevLanguage::Parser::executeBaseVariable(void)
{
    if (base_variable_expr != NULL)
    {
        base_variable = base_variable_expr->evaluateContent(Workspace::userWorkspace());
    }
}

/** 
 * Give flex a line to parse
 */
void RevLanguage::Parser::getline(char* buf, size_t maxsize)
{

    if (!rrcommand.good())
    {
        foundEOF = true;
        //        buf[0] = EOF;
        //        buf[1] = '\0';
        buf[0] = '\0';
    }
    else
    {
        foundNewline = false;
        rrcommand.getline(buf, long(maxsize) - 3);
        // Deal with line endings in case getline uses non-Unix endings
        size_t i = strlen(buf);
        if (i >= 1 && buf[i - 1] == '\r')
        {
            buf[i - 1] = '\n';
        }
        else if (i >= 2 && buf[i - 1] == '\n' && buf[i - 2] == '\r')
        {
            buf[i - 2] = '\n';
            i--;
        }
        else if (i == 0 || (i >= 1 && buf[i - 1] != '\n'))
        {
            buf[i++] = '\n';
        }
        buf[i] = '\0';
        
    }

}

/** This function gets help info about a symbol */
int RevLanguage::Parser::help(const std::string& symbol) const
{
    
    // Get some help
    RevBayesCore::RbHelpSystem& hs = RevBayesCore::RbHelpSystem::getHelpSystem();
    if ( hs.isHelpAvailableForQuery(symbol) )
    {
        const RevBayesCore::RbHelpEntry& h = hs.getHelp( symbol );
        RevBayesCore::HelpRenderer hRenderer;
        std::string hStr = hRenderer.renderHelp(h, RbSettings::userSettings().getLineWidth() - RevBayesCore::RbUtils::PAD.size());
        UserInterface::userInterface().output("\n", true);
        UserInterface::userInterface().output("\n", true);
        UserInterface::userInterface().output(hStr, true);
    }
    else
    {
        
        RBOUT("Help is not available for \"" + symbol + "\"");
        
    }
    
    // Return success
    return 0;
}


/** This function gets help info about a symbol */
int RevLanguage::Parser::help(const std::string& baseSymbol, const std::string& symbol) const
{
    
    // Get some help
    RevBayesCore::RbHelpSystem& hs = RevBayesCore::RbHelpSystem::getHelpSystem();
    if ( hs.isHelpAvailableForQuery( baseSymbol, symbol) )
    {
        const RevBayesCore::RbHelpEntry& h = hs.getHelp( baseSymbol, symbol );
        RevBayesCore::HelpRenderer hRenderer;
        std::string hStr = hRenderer.renderHelp(h, RbSettings::userSettings().getLineWidth() - RevBayesCore::RbUtils::PAD.size());
        UserInterface::userInterface().output("\n", true);
        UserInterface::userInterface().output("\n", true);
        UserInterface::userInterface().output(hStr, true);
    }
    else
    {
        
        RBOUT("Help is not available for \"" + baseSymbol + "." + symbol + "\"");
        
    }
    
    // Return success
    return 0;
}


/**
 * This function prints help info about a function if it sees a function call.
 * The function is called from the bison code, which is responsible for
 * deleting the syntax tree (the function call).
 */
int RevLanguage::Parser::help(const SyntaxFunctionCall* root) const
{
    
    RlString symbol = root->getFunctionName();
    
    return help( symbol.getValue() );
}


extern int yyparse(void); // Defined in grammar.tab.cpp (from gammar.y)

/**
 * Set parser mode. Use this function to set the parser mode. Use
 * Parser::CHECKING for state checking and Parser::PROCESSING for
 * command processing.
 */
void RevLanguage::Parser::setParserMode(ParserMode mode)
{
    base_variable = NULL;
    function_name = "";
    argument_label = "";
    base_variable_expr = NULL;

    parser_mode = mode;
}

/**
 * Process command buffer with the help of the bison-generated code.
 *
 * @param   command     command RlString; modified to contain any unprocessed incomplete statement residue, "" if none
 * @return  integer flag indicating status:
 *          1 - (last) statement in command buffer is incomplete
 *          2 - syntax error encountered
 *          0 - all statements in command buffer syntactically correct, complete and executed
 *
 * @note If an exception occurs on execution of a syntactically correct
 *       statement, the exception msg is printed here and the return
 *       signal is set to 2. Any remaining part of the command buffer
 *       is discarded.
 */
int RevLanguage::Parser::processCommand(std::string& command, Environment* env)
{

    // make sure mode is not checking
    this->setParserMode(EXECUTING);

    executionEnvironment = env;

    // Break command into Rev lines
    std::vector<std::string> lines;
    try
    {
        breakIntoLines(command, lines, true);
    }
    catch (RbException& rbException)
    {

        // Print message
        std::ostringstream msg;
        rbException.print(msg);
        msg << std::endl;
        RBOUT(msg.str());

        // We printed a message so we dealt with the problem
        // We exit immediately, discarding any command buffer content
        // We return 2 to signal a problem, which the caller may choose to ignore or act upon
        command = "";
        return 2;
    }

    // Initialize flex column count
    yycolumn = 1;

    // Call Bison code, which calls Flex code, which calls rrinput
    for (std::vector<std::string>::iterator i = lines.begin(); i != lines.end(); i++)
    {

        /* prepare globals for call to parser */
        rrcommand.str((*i));
        rrcommand.clear();
        foundNewline = false;
        foundEOF = false;
        foundErrorBeforeEnd = false;

#if 0
        // NOTE! This code is only for testing the CHECKING mode of the parser
        // DO NOT UNCOMMENT IN NORMAL USE
        /* Set parser mode into checking */
        setParserMode(CHECKING);
#endif


        int result;
        try
        {
            result = yyparse();
        }
        catch (RbException& rbException)
        {

            // Catch a quit request in case it was not caught before
            if (rbException.getExceptionType() == RbException::QUIT)
            {
                Workspace::userWorkspace().clear();
                Workspace::globalWorkspace().clear();
                
#ifdef RB_MPI
                MPI_Finalize();
#endif
                
                exit(0);
            }
            // All other uncaught exceptions
            
            rbException.print(std::cerr);
            std::cerr << std::endl;
            
            // We exit immediately, discarding any remaining buffer content
            // We return 2 to signal a problem, which the caller may choose to ignore or act upon
            command = "";
            return 2;
        }

#if 0
        // NOTE! This code is only for testing the CHECKING mode of the parser
        // DO NOT UNCOMMENT IN NORMAL USE
        if (result == 1) {
            std::cerr << "Function name is: " << function_name << std::endl;
            std::cerr << "Argument label is: " << argument_label << std::endl;
            if (base_variable == NULL)
                std::cerr << "No base variable" << std::endl;
            else
                std::cerr << "Base variable is: " << base_variable->getName() << std::endl;
        }
#endif

        if (result == 0)
        {

        }
        else if (result == 2)
        {

            // The execute function printed the error message resulting from the exception
            // We exit immediately, discarding any remaining buffer content
            // We return 2 to signal a problem, which the caller may choose to ignore or act upon
            command = "";
            return 2;
        }
        else if (foundErrorBeforeEnd == true)
        {

            std::ostringstream msg;
            if (yylloc.first_column == yylloc.last_column)
            {
                msg << "Syntax error while reading character '" << command[size_t(yylloc.first_column - 1)] << "' at position " << yylloc.first_column << " in command:" << std::endl;
            }
            else
            {
                msg << "Syntax error while reading characters \"";
                for (int j = yylloc.first_column; j <= yylloc.last_column; ++j)
                {
                    msg << command[size_t(j - 1)];
                }
                msg << "\" at position " << yylloc.first_column << " to " << yylloc.last_column << " in command:" << std::endl;
            
            }
            msg << command;

            RBOUT(msg.str());
            command = "";
            return 2;
        }
        else if (foundNewline == true && foundEOF == false)
        {

            /* Replace newline with space */
            std::string temp = (*i);
            temp[temp.size() - 1] = ' ';

            std::vector<std::string>::iterator j = i;
            j++;
            if (j == lines.end())
            {

                /* If no more input lines, we need to ask for more */
                command = temp;
                return 1;
            }
            else
            {

                /* If more input lines, put temp in front before proceeding */
                (*j) = temp + (*j);
            }
        }
        else if (foundNewline == true && foundEOF == true)
        {

            std::vector<std::string>::iterator j = i;
            j++;
            if (j == lines.end())
            {

                /* If no more input lines, we need to ask for more */
                command = (*i);
                return 1;
            } else {

                /* If more input lines, put current line in front before proceeding */
                (*j) = (*i) + (*j);
            }
        }
        else
        {

            RBOUT("Unknown parse error");
            command = "";
            return 2;
        }
        
    }

    /* Successfully processed all statements in command buffer */
    command = "";
    return 0;
}

ParserInfo Parser::checkCommand(std::string& command, Environment* env)
{

    setParserMode(CHECKING);
    executionEnvironment = env;
    ParserInfo pi;
    pi.message = "";

    // Break command into Rev lines
    std::vector<std::string> lines;
    try
    {
        pi = breakIntoLines(command, lines, false);
    }
    catch (RbException& rbException)
    {
        pi.message.append("Caught an exception while breaking command buffer into lines\n\r");
        // unrecoverable error        
        // We return 2 to signal a problem, which the caller may choose to ignore or act upon
        pi.result = 2;
        return pi;
    }

    // Initialize flex column count
    yycolumn = 1;

    // Call Bison code, which calls Flex code, which calls rrinput
    for (std::vector<std::string>::iterator i = lines.begin(); i != lines.end(); i++) {
        pi.function_name = function_name;
        pi.base_variable = NULL;
        pi.argument_label = argument_label;
        pi.result = 0;


        /* prepare globals for call to parser */
        rrcommand.str((*i));
        rrcommand.clear();
        foundNewline = false;
        foundEOF = false;
        foundErrorBeforeEnd = false;

        pi.message.append("\n\rCalling bison with rrcommand:\n\r").append(rrcommand.str()).append("\n\r");

        int result;
        try {
            result = yyparse();
        }        catch (RbException& rbException) {
            pi.message.append("Caught an exception calling yyparse\n\r");
            // Catch a quit request in case it was not caught before
            if (rbException.getExceptionType() == RbException::QUIT) {
                pi.message.append("ignoring QUIT request");
            }
            // All other uncaught exceptions
            pi.message.append("Abnormal exception during parsing or execution of statement; discarding any remaining command buffer\n\r");

            // try parsing next line             
            pi.result = 2;
            continue;
        }

        if (result == 0 || result == 1) {

            pi.function_name = function_name;
            pi.base_variable = base_variable;
            pi.argument_label = argument_label;

            // valid result, nothing more to do here?
            //break;

            // I think we should parse all lines
            continue;
        }

        if (foundNewline == true && foundEOF == false) {
            pi.message.append("Incomplete statement ending with inappropriate newline; replaced with space\n\r");

            /* Replace newline with space */
            std::string temp = (*i);
            temp[temp.size() - 1] = ' ';

            std::vector<std::string>::iterator j = i;
            j++;
            if (j == lines.end()) {
                /* If no more input lines, we need to ask for more */
                pi.message.append("Reached end; asking for more content to append to (last) statement.\n\r\n\r");
                command = temp;
                pi.result = 1;
                return pi;
            } else {
                /* If more input lines, put temp in front before proceeding */
                pi.message.append("Continuing with next Rev line in command buffer.\n\r\n\r");
                (*j) = temp + (*j);
            }
        } else if (foundNewline == true && foundEOF == true) {
            pi.message.append("Incomplete statement ending with appropriate newline.\n\r");
            std::vector<std::string>::iterator j = i;
            j++;
            if (j == lines.end()) {
                /* If no more input lines, we need to ask for more */
                pi.message.append("Reached end; asking for more content to append to (last) statement.\n\r\n\r");
                command = (*i);
                pi.result = 1;
                return pi;
            } else {
                /* If more input lines, put current line in front before proceeding */
                pi.message.append("Continuing with next Rev line in command buffer.\n\r\n\r");
                (*j) = (*i) + (*j);
            }
        }
    }

    // make sure executing mode is set by default
    setParserMode(EXECUTING);
    return pi;
}

/** Global call-back function for the flex-generated code in lex.yy.cpp */
void rrinput(char* buf, int& result, size_t maxsize) {

    buf[0] = '\0';
    RevLanguage::Parser::getParser().getline(buf, maxsize);
    result = strlen(buf);
}
