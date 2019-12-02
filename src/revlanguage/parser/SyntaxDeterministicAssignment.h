#ifndef SyntaxDeterministicAssignment_H
#define SyntaxDeterministicAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Deterministic assignments
     *
     * This syntax element class handles deterministic assignments,
     * sometimes called equation assignments. These are assignment
     * statements of the type:
     *
     *    a := foo( x )
     *
     * using the colon-equal or equation assignment operator.
     *
     * The intention of deterministic assignments is to create
     * deterministic variables for inclusion in DAGs. These are
     * typically used for variable conversions, but also for
     * operations involving the construction of complex objects
     * from component parts.
     *
     * Unlike constant variables, deterministic variables change
     * their value when variables included in the right-hand side
     * expression change their value or are reassigned new values.
     *
     * For this reason, the right-hand side expression of a deterministic
     * assignment is referred to as a dynamic expression. Its semantic
     * value is evaluated using the evaluateContent(dynamic=true) function,
     * instead of the standard evaluateContent(dynamic=false) function.
     *
     * Many syntax elements evaluate to the same objects in both contexts,
     * for instance constants. For other syntax elements, the behavior
     * changes subtly but with important consequences. This is most
     * obvious in the variable syntax element but also in the function
     * call, unary and binary expression elements.
     */
    class SyntaxDeterministicAssignment : public SyntaxAssignment {
        
    public:
        SyntaxDeterministicAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                                      //!< Basic constructor

	    virtual                             ~SyntaxDeterministicAssignment();                                               //!< Destructor
        
        // Basic utility functions
        SyntaxDeterministicAssignment*      clone() const;                                                                  //!< Clone object
                
    protected:
        
        void                                assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);                           //!< The assignment operation.
        bool                                isDynamic(void);                                                                //!< Should the rhs be evaluated dynamically?
    
    };
    
}

#endif

