#ifndef SyntaxBinaryExpr_H
#define SyntaxBinaryExpr_H

#include "SyntaxElement.h"

#include <iostream>
#include <vector>


namespace RevLanguage {

    /**
     * @brief Binary expression syntax element
     *
     * This syntax element is used to hold binary expressions. Member
     * variables are the left-hand side and right-hand side operands,
     * and the operator type.
     *
     * Binary expressions are expressions of the type
     *
     *    a + b
     *
     * and they effectively correspond to calls of internal bivariate
     * functions.
     *
     * Just like function calls, binary expressions evaluate in
     * differently ways in dynamic and static expression contexts. For
     * more explanation of this, see the function call syntax element.
     *
     * The operand arguments of binary expressions are understood to
     * be dynamic arguments. See the function call syntax element
     * for more details on different argument types.
     */
    class SyntaxBinaryExpr : public SyntaxElement {

    public:
        // Binary operator types
        enum                        operatorT { Range, Add, Sub, Mul, Div, Mod, Exp, Lt, Le, Eq, Ne, Ge, Gt, And, Or, And2, Or2 };

        static std::string          opCode[];                                                                   //!< Operator codes for printing

        // Constructors and destructor
        SyntaxBinaryExpr(SyntaxBinaryExpr::operatorT op, SyntaxElement* lhs, SyntaxElement* rhs);               //!< Standard constructor
        SyntaxBinaryExpr(const SyntaxBinaryExpr& x);                                                            //!< Copy constructor

        virtual                    ~SyntaxBinaryExpr();                                                         //!< Destroy operands

        // Assignment operator
        SyntaxBinaryExpr&           operator=(const SyntaxBinaryExpr& x);                                       //!< Assignment operator

        // Basic utility functions
        SyntaxElement*              clone() const;                                                              //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);  //!< Get semantic value
        bool                        isConstExpression(void) const;                                              //!< Is the expression constant?
        bool                        isFunctionSafe(const Environment&       env,
                                                   std::set<std::string>&   localVars) const;                   //!< Is this element safe in a function?

    protected:
        SyntaxElement*              left_operand;                                                               //!< The left operand
        SyntaxElement*              right_operand;                                                              //!< The right operand
        enum operatorT              operation;                                                                  //!< The type of operation
    
    };
    
}

#endif

