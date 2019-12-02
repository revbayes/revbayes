#ifndef SyntaxConstantAssignment_H
#define SyntaxConstantAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Constant assignment
     *
     * Objects of this class are used to hold left-arrow or constant
     * assignment statements. The member variables are the lhs
     * (left-hand side) and rhs (righ-hand side) expressions.
     *
     * A constant assignment or left-arrow assignment is a statement
     * of the type
     *
     *    a <- x + y
     *
     * and uses a left-arrow assignment operator. The rhs expression
     * is evaluated as a constant or static expression using the
     * standard evaluateContent() function, and not its dynamic or
     * lhs versions.
     */
    class SyntaxConstantAssignment : public SyntaxAssignment {
        
    public:
        SyntaxConstantAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                                   //!< Basic constructor
        
        virtual                    ~SyntaxConstantAssignment();                                                     //!< Destructor
        
        
        // Basic utility functions
        SyntaxConstantAssignment*   clone() const;                                                                  //!< Clone object
        
    protected:
        
        void                        assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);                           //!< The assignment operation.
        
    };
    
}

#endif

