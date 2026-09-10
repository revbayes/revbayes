#ifndef SyntaxAssignment_H
#define SyntaxAssignment_H

#include "SyntaxElement.h"

#include <iostream>

namespace RevLanguage {
    
    /**
     * @brief Syntax assignment: base class for all assignment types (constant, deterministc and stochastic)
     *
     * This is an abstract base class for all assignment operations. Here we provide some common functionality
     * that all assignment share.
     */
    class SyntaxAssignment : public SyntaxElement {
        
    public:
        SyntaxAssignment(SyntaxElement* lhsExpr, SyntaxElement* rhsExpr);                                           //!< Basic constructor
        SyntaxAssignment(const SyntaxAssignment& x);                                                                //!< Copy constructor
        
	    virtual                    ~SyntaxAssignment();                                                             //!< Destructor
        
        // Assignment operator
        SyntaxAssignment&           operator=(const SyntaxAssignment& x);                                           //!< Assignment operator
        
        // Basic utility functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value
        bool                        isAssignment(void) const;                                                       //!< Is this an assignment statement?
        bool                        isFunctionSafe(const Environment&       env,
                                                   std::set<std::string>&   localVars) const;                       //!< Is this element safe in a function?
        
        // pure virtual regular functions
        virtual SyntaxAssignment*   clone() const = 0;                                                              //!< Clone object
        
    protected:
        // protected methods
        virtual void                assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs) = 0;                       //!< The assignment operation.
        virtual bool                isDynamic(void);                                                                //!< Should the rhs be evaluated dynamically?
//        void                        removeElementVariables(Environment &env, RevPtr<RevVariable> &the_var);             //!< Removing element variables from this variable.
        
        // protected members
        SyntaxElement*              lhsExpression;                                                                  //!< The lhs expression
        SyntaxElement*              rhsExpression;                                                                  //!< The rhs expression
        
    };
    
}

#endif

