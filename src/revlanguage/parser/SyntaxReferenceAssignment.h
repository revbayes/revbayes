#ifndef SyntaxReferenceAssignment_H
#define SyntaxReferenceAssignment_H

#include "SyntaxElement.h"
#include "SyntaxVariable.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Reference variable assignment
     *
     * This syntax element class is used to hold reference variable
     * assignment statements. The member variables are the lhs
     * (left-hand side) and rhs (righ-hand side) expressions.
     *
     * A reference variable assignment is a statement like
     *
     *    a <-& b
     *
     * using the reference assignment operator ('<-&'). The assignment
     * results in 'a' being an alias for 'b'. That is, 'a' will be bound
     * to the variable 'b', so that it always has the same value as 'b'.
     * It will act as a synonym for 'b' in all contexts except in
     * basic assignments (constant, control, stochastic, and deterministic
     * assignments).
     * 
     * Unlike some languages, where variable aliases are permanent
     * once they have been created, the Rev reference variables can
     * be redefined through assignment. That is, if you issed the
     * statements
     *
     *    b <- 1
     *    a <-& b
     *    a <- 2
     *
     * then a would initially have the value 1, the same as b, but
     * after the final assignment, a will have the value 2 and b
     * will have the value 1. This behavior is more in alignment with
     * the replacement assignment philosophy of the rest of the Rev
     * language.
     *
     * With a reference variable you can call all operations that
     * change the value of b in place. That is, in this series of
     * statements
     *
     *    b ~ norm( 0, 1 )
     *    a <-& b
     *    a.setValue( 2.1 )
     *
     * the final statement will result in the value of 'b' (and thus
     * of 'a') being set to 2.1.
     *
     * Similarly, the value of the referenced variable can be modified
     * with increment and decrement operators, as well as with addition,
     * subtraction, multiplication and division assignments. Thus,
     *
     *    a <-& b
     *    a++;
     *
     * will result in 'b' being increased by 1.
     *
     * Reference variables are important in some contexts, particularly
     * in user-defined moves, where the move parameter needs to be
     * a reference variable to the random variable being affected by
     * the move.
     */
    class SyntaxReferenceAssignment : public SyntaxElement {
        
    public:
        SyntaxReferenceAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                      //!< Basic constructor
        SyntaxReferenceAssignment(const SyntaxReferenceAssignment& x);                                  //!< Copy constructor

	    virtual                         ~SyntaxReferenceAssignment();                                   //!< Destructor
        
        // Assignment operator
        SyntaxReferenceAssignment&      operator=(const SyntaxReferenceAssignment& x);                  //!< Assignment operator
        
        // Basic utility functions
        SyntaxReferenceAssignment*      clone() const;                                                  //!< Clone object
        bool                            isAssignment(void) const;                                       //!< Is this an assignment statement?
        
        // Regular functions
        RevPtr<RevVariable>             evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);          //!< Get semantic value
        bool                            isFunctionSafe(const Environment&       env,
                                                       std::set<std::string>&   localVars) const;       //!< Is this element safe in a function?
        
    protected:
        
        SyntaxElement*                  lhsExpression;                                                  //!< The lhs expression
        SyntaxElement*                  rhsExpression;                                                  //!< The rhs expression
        
    };
    
}

#endif

