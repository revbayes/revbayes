#ifndef SyntaxStatement_H
#define SyntaxStatement_H

#include "SyntaxElement.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    /**
     * @brief Statement syntax element
     *
     * This syntax element is used to store and execute statements
     * like control structures, for loops, if-else statements etc.
     */
    class SyntaxStatement : public SyntaxElement {

    public:
        // Statement types
        enum                        statementT { If, IfElse, For, While, Next, Break, Return };                         //!< Statement type codes
        static std::string          stmtName[];                                                                         //!< Statement names for printing

        // Constructors and destructor
                                                    SyntaxStatement(statementT type);                                               //!< NEXT, BREAK, RETURN statements
                                                    SyntaxStatement(SyntaxStatement::statementT type,
                                                                    SyntaxElement*              expr);                              //!< RETURN statement
                                                    SyntaxStatement(statementT                  type,
                                                                    SyntaxElement*              cond,
                                                                    std::list<SyntaxElement*>*  stmts);                             //!< IF, FOR, WHILE statements
                                                    SyntaxStatement(statementT                  type,
                                                                    SyntaxElement*              expr,
                                                                    std::list<SyntaxElement*>*  stmts1,
                                                                    std::list<SyntaxElement*>*  stmts2);                            //!< IF_ELSE statements
                                                    SyntaxStatement(const SyntaxStatement& x);                                      //!< Copy constructor
	    virtual                                    ~SyntaxStatement();                                                              //!< Destroy operands

        // Assignment operator
        SyntaxStatement&                            operator=(const SyntaxStatement& x);                                            //!< Assignment operator

        // Basic utility functions
        SyntaxStatement*                            clone() const;                                                                  //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>                         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value

    protected:
        bool                                        isTrue(SyntaxElement* expression, const std::shared_ptr<Environment>& env) const;  //!< Does expression evaluate to true?

        enum statementT                             statementType;                                                                  //!< The type of statement
        SyntaxElement*                              expression;                                                                     //!< Expression, conditional expr, or for condition
        std::list<SyntaxElement*>*                  statements1;                                                                    //!< First set of statements
        std::list<SyntaxElement*>*                  statements2;                                                                    //!< Second set of statements
    
    };

}

#endif

