#ifndef SyntaxStochasticAssignment_H
#define SyntaxStochasticAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Stochastic assignment
     *
     * Objects of this class are used to hold tilde or stochastic
     * assignment statements. The member variables are the lhs
     * (left-hand side) and rhs (righ-hand side) expressions.
     *
     * A stochastic assignment or tilde assignment is a statement
     * of the type
     *
     *    a ~ norm( 0.0, 1.0 )
     *
     * and uses a tilde assignment operator.
     *
     * The rhs expression is evaluated as a dynamic expression.
     * The resulting variable is checked for its type, to make sure
     * that it is a distribution. A random variable is then drawn
     * from the distribution, and it is inserted into the variable
     * resulting from evaluation of the left-hand side expression.
     */
    class SyntaxStochasticAssignment : public SyntaxAssignment {
        
    public:
        SyntaxStochasticAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                                 //!< Basic constructor
        
	    virtual                        ~SyntaxStochasticAssignment();                                               //!< Destructor
        
        // Basic utility functions
        SyntaxStochasticAssignment*     clone() const;                                                              //!< Clone object
        
    protected:
        
        void                            assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);                       //!< The assignment operation.
        bool                            isDynamic(void);                                                            //!< Should the rhs be evaluated dynamically?

    };
    
}

#endif

