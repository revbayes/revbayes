#ifndef SyntaxWorkspaceVariableAssignment_H
#define SyntaxWorkspaceVariableAssignment_H

#include "SyntaxAssignment.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Workspace variable assignment
     *
     * This syntax element class is used to hold workspace variable
     * assignment statements. The member variables are the lhs
     * (left-hand side) and rhs (righ-hand side) expressions.
     *
     * A workspace variable assignment is a statement like
     *
     *    i = 1
     *    m = mvScale( xxx )
     *
     * using the equal assignment operator.
     *
     * A workspace variable assignment differs from a constant
     * assignment in two ways. First, workspace objects (moves, monitors, etc.)
     * can only be created using a workspace variable assignment. 
     * Second, model variables (integer, etc.) created using the workspace
     * variable assignment will not be included in any DAGs even if it occurs
     * in a dynamic (deterministic or stochastic) expression.
     *
     * The workspace variable behavior is typically used for loop
     * variables. The for loop variable is implicitly set to be
     * a workspace variable.
     */
    class SyntaxWorkspaceVariableAssignment : public SyntaxAssignment {
        
    public:
        SyntaxWorkspaceVariableAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                          //!< Basic constructor
        
        virtual                                 ~SyntaxWorkspaceVariableAssignment();                               //!< Destructor
        
        // Basic utility functions
        SyntaxWorkspaceVariableAssignment*      clone() const;                                                      //!< Clone object
        
    protected:
        
        void                                    assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs);               //!< The assignment operation.
        bool                                    isDynamic(void);                                                            //!< Should the rhs be evaluated dynamically?

    };
    
}

#endif

