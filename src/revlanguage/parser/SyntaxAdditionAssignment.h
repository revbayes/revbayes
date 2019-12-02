#ifndef SyntaxAdditionAssignment_H
#define SyntaxAdditionAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Addition assignment operator ('a += b')
     *
     * The addition assignment operator adds the value of the right-hand side variable
     * to the value of the left-hand side variable.
     *
     * Note that the left-hand side variable must be either a constant variable or a 
     * control variable, and that the addition assignment is interpreted as a constant
     * assignment or a control variable assignment depending on the original type of the
     * left-hand side variable. Thus, if the left-hand side variable is a constant, the
     * addition assignment statement 'a += b' is equivalent to the statement:
     *
     *    a <- a + b
     *
     * while if the left-hand side variable is a control variable, the addition assignment
     * is equivalent to
     *
     *    a <<- a + b
     *
     * Note that an addition assignment could not be a deterministic assignment, because
     * then it would be rejected as creating loops in the DAG, i.e., 'a' would be a parent
     * of itself.
     */
    class SyntaxAdditionAssignment : public SyntaxAssignment {
        
    public:
        SyntaxAdditionAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                           //!< Standard constructor
	    virtual                             ~SyntaxAdditionAssignment();                                    //!< Destructor
        
        // Basic utility functions
        SyntaxAdditionAssignment*           clone() const;                                                  //!< Clone object
        
    protected:
        void                                assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);           //!< The assignment operation.
        
    };
    
}

#endif

