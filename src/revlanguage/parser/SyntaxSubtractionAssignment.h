#ifndef SyntaxSubtractionAssignment_H
#define SyntaxSubtractionAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Subtraction assignment operator ('a -= b')
     *
     * The subtraction assignment operator subtracts the value of the right-hand side variable
     * from the value of the left-hand side variable.
     *
     * Note that the left-hand side variable must be either a constant variable or a
     * control variable, and that the subtraction assignment is interpreted as a constant
     * assignment or a control variable assignment depending on the original type of the
     * left-hand side variable. Thus, if the left-hand side variable is a constant, the
     * subtraction assignment statement 'a -= b' is equivalent to the statement:
     *
     *    a <- a - b
     *
     * while if the left-hand side variable is a control variable, the subtraction assignment
     * is equivalent to
     *
     *    a <<- a - b
     *
     * Note that a subtraction assignment could not be a deterministic assignment, because
     * then it would be rejected as creating loops in the DAG, i.e., 'a' would be a parent
     * of itself.
     */
    class SyntaxSubtractionAssignment : public SyntaxAssignment {
        
    public:
        SyntaxSubtractionAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                        //!< Standard constructor
	    virtual                            ~SyntaxSubtractionAssignment();                                  //!< Destructor
        
        // Basic utility functions
        SyntaxSubtractionAssignment*        clone() const;                                                  //!< Clone object
        
    protected:
        
        void                                assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);           //!< The assignment operation.
        
    };
    
}

#endif

