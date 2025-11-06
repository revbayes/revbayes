#ifndef Parser_H
#define Parser_H

#include <cstddef>
#include <list>
#include <string>
#include <sstream>
#include <vector>

#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {


    class ParserInfo {
    public:
        bool inQuote;
        bool inComment;
        std::vector<std::string> lines;
        int result;
        std::string message;
        std::string function_name;
        std::string argument_label;
        RevPtr<RevVariable> base_variable;
    };
    
    class Environment;
    class SyntaxElement;
    class SyntaxFunctionCall;
    class SyntaxVariable;

    /**
     * @brief Singleton Parser class
     *
     * The parser relies on automatically generated code produced
     * by bison and flex, in the files grammar.tab.cpp, grammar.tab.h
     * and lex.yy.cpp.
     *
     * You use the parser by calling the processCommand function with
     * a string buffer containing one or more, possibly incomplete, Rev
     * lines. The function breaks the buffer into Rev lines using the
     * breakIntoLines help function. Then it cycles through the lines,
     * using the bison and flex code through the call to bison's yyparse
     * function. The call to yyparse generates one or more calls to flex's
     * yylex function, which in turn calls the global rrinput function,
     * defined in the Parser.cpp file. The rrinput function uses the parser's
     * getline function to get a single line of input from the global
     * stringstream rrcommand, which is loaded with one Rev line at a time by
     * the processCommand function before the call to yyparse.
     * 
     * The call to yyparse generates an error code if there is a syntax error
     * in the command buffer. The destructors in the bison-generated code should
     * deal with the deletion of the syntax element subtrees that have been created
     * up to the point when the error was encountered.
     * 
     * For each complete and syntactically correct statement found in the buffer, the
     * bison code generates a call to the execute or help functin of the parser once
     * the syntax element tree is complete. The execute function executes the syntax
     * tree corresponding to the statement, while the help function gets appropriate help
     * information if the tree corresponds to a single statement that is a help request.
     * This behavior is defined in the grammar.y file. Once the syntax tree has been
     * evaluated, it is deleted by the execute or the help function in the parser.
     *
     * Note that there can be nested calls to the parser, for instance
     * if the 'source' function starts parsing a file as part of the
     * execution of a statement resulting from the parsing of a
     * command buffer containing a 'source' call. Rather than using
     * multiple parsers to handle this, each with its own flex buffer, we
     * use a unique global parser. This works because we rely on the fact that
     * each call to processCommand generates its own command buffer with its
     * constituent Rev lines. By making sure we only call the bison and flex
     * code with one Rev line at a time, we avoid any issues with global state
     * flags or buffers that have to be maintained across calls to yyparse.
     * Nested calls to processCommand work because each call maintains its own set
     * of Rev lines. Once an inner call to processCommand has finished, the outer
     * call can continue processing its remaining lines, if any.
     *
     * If multiple Rev lines had been fed to a single flex buffer, then the remaining
     * lines of the outer call would have been discarded by the time the inner call had
     * finished. For instance, this problem would occur with a command string
     * such as 'source("somefile"); <this statement would be discarded>'.
     *
     * The parser class also supports processing of a partial command line, to find out
     * the parser state at the end of the command line. Local variables, such as
     * function_name, base_variable_expression etc, as well as some of the public functions
     * support this use case.
     */
        
    class Parser {

    public:
        
        enum ParserMode { CHECKING, EXECUTING };                                                    //!< Parser modes

        // Regular functions
        int                 execute(SyntaxElement* root, const std::shared_ptr<Environment>& env) const;  //!< Execute the syntax tree
        void                getline(char* buf, size_t maxsize);                                     //!< Give flex one line to process
        int                 help(const std::string& symbol) const;                                  //!< Get help for a symbol
        int                 help(const std::string& baseSymbol, const std::string& symbol) const;   //!< Get help for a symbol
        int                 help(const SyntaxFunctionCall* root) const;                             //!< Get help for a function call
        int                 processCommand(std::string& command, const std::shared_ptr<Environment>& env);     //!< Process command with help from Bison

        // State checking functions
        bool                isChecking(void) { return parser_mode == CHECKING; }                     //!< Are we in state-checking mode?
        void                setBaseVariable(SyntaxVariable* var) { base_variable_expr = var; }        //!< Set base variable expression
        void                executeBaseVariable(void);                                              //!< Execute base variable expression
        void                setfunction_name( const std::string& n ) { function_name = n; }           //!< Set function name
        void                setArgumentLabel( const std::string& n ) { argument_label = n; }         //!< Set argument label
        
        ParserInfo          checkCommand(std::string& command, const std::shared_ptr<Environment>& env);     //!< Parse command without executing it
        
        /** Get singleton parser */
        static Parser& getParser( void ) {
            static Parser theParser;
            return theParser;
        }
        
   private:
                                Parser(void);                                                           //!< Prevent public construction
                                Parser(const Parser& x) {}                                              //!< Prevent copy construction

        Parser&                 operator=(const Parser& w) { return (*this); }                          //! Prevent assignment

        // Help functions
        ParserInfo              breakIntoLines(const std::string& cmd, std::vector<std::string>& lines, bool validate) const;    //!< Break a command string buffer into Rev lines
        ParserInfo              breakIntoLines(const std::string& cmd, std::vector<std::string>& lines) const;    //!< Break a command string buffer into Rev lines
        void                    setParserMode(ParserMode mode);                                         //!< Set the parser mode
        
        // Member variables
        ParserMode              parser_mode;                                             //!< The current parser mode
        std::string             function_name;                                           //!< Function name of end state, if any
        std::string             argument_label;                                          //!< Argument label of end state, if any
        SyntaxVariable*         base_variable_expr;                                       //!< Base variable expression of end state
        RevPtr<RevVariable>     base_variable;                                           //!< Base variable of end state

    };

}

extern std::shared_ptr<RevLanguage::Environment> executionEnvironment;

// Global call-back function for flex-generated code
void rrinput(char* buf, int& result, size_t maxsize);

// Global command buffer string stream for rrinput
extern std::stringstream rrcommand;

// Global flags indicating flex state
extern bool foundNewline;
extern bool foundErrorBeforeEnd;
extern bool foundEOF;

// Global column and token length variables for flex
extern size_t yycolumn;
extern int    yyleng;

#endif
