#ifndef SyntaxDivisionAssignment_H
#define SyntaxDivisionAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Division assignment operator ('a /= b')
     *
     * The division assignment operator divides the value the left-hand side variable
     * with the right-hand side variable.
     *
     * Note that the left-hand side variable must be either a constant variable or a
     * control variable, and that the division assignment is interpreted as a constant
     * assignment or a control variable assignment depending on the original type of the
     * left-hand side variable. Thus, if the left-hand side variable is a constant, the
     * division assignment statement 'a /= b' is equivalent to the statement:
     *
     *    a <- a / b
     *
     * while if the left-hand side variable is a control variable, the addition assignment
     * is equivalent to
     *
     *    a <<- a / b
     *
     * Note that a division assignment could not be a deterministic assignment, because
     * then it would be rejected as creating loops in the DAG, i.e., 'a' would be a parent
     * of itself.
     */
    class SyntaxDivisionAssignment : public SyntaxAssignment {
        
    public:
        SyntaxDivisionAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                           //!< Standard constructor
	    virtual                             ~SyntaxDivisionAssignment();                                    //!< Destructor
        
        // Basic utility functions
        SyntaxDivisionAssignment*           clone() const;                                                  //!< Clone object

    protected:
        
        void                                assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);           //!< The assignment operation.
        
    };
    
}

#endif

