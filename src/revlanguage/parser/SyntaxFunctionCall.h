#ifndef SyntaxFunctionCall_H
#define SyntaxFunctionCall_H

#include "SyntaxElement.h"
#include "SyntaxLabeledExpr.h"
#include "SyntaxVariable.h"

#include <iostream>
#include <list>
#include <utility>


namespace RevLanguage {

    /**
     * @brief Function call syntax element
     *
     * This syntax element is used to hold function calls. We store the arguments,
     * function name, and (if applicable) the variable of which the function is a member.
     * 
     * Function calls are expressions of the type
     *
     *    foo( x )
     *
     * or
     *
     *    a.foo( x )
     *
     * as might be expected. Note that expressions like 'f(x)[1]' are parsed out as
     * variable elements and not function calls. In such cases, however, the variable
     * syntax element will contain a function call syntax element holding the call 'f(x)'.
     * If the function is a member function, as in 'a.foo( x )', then 'a' is referred to
     * as the base variable, stored in base_variable. If there is no base variable, 
     * base_variable is set to NULL.
     *
     * The function call will (potentially) behave differently depending on whether
     * or not it occurs in a dynamic expression. Recall that dynamic expressions are
     * stochastic and deterministic expression on the right-hand side of stochastic
     * and deterministic assignments.
     *
     * If a function call occurs in a static (constant) expression, it will be evaluated
     * once, and the return type is a constant or composite variable. If the call occurs
     * in a dynamic expression, then the call results in the creation of a deterministic
     * variable of the return type, the value of which will change dynamically if the
     * arguments passed in are changed, or if one reassigns to the variables in those
     * expressions.
     *
     * The dynamic and static evaluation of function calls are performed using the 
     * functions evaluateDynamicContent() and evaluateContent(), respectively.
     *
     * Function arguments can be specified to be 'const', 'dynamic' or 'mutable'. By default,
     * arguments are 'dynamic', which means that they are evaluated as dynamic expressions
     * if the function call itself is part of a dynamic expression. A 'const' argument is
     * always evaluated as a constant (static) expression regardless of whether the function
     * call occurs in a dynamic or a static expression.
     *
     * Neither 'const' nor 'dynamic' arguments can be modified by the function. A 'mutable'
     * argument, however, is evaluated as a dynamic expression and the resulting variable
     * can be modified by the function. Functions taking 'mutable' arguments cannot be used
     * in DAGs because they potentially create hard-to-detect loops in the DAG structure
     * if they use their power to change the arguments passed in.
     *
     * The argument matching rules in Rev are similar to those in R, but the fact that Rev
     * is a typed language presents some additional complexity. See RevLanguage::Function
     * for more detailed explanation of the argument matching in Rev.
     */
    class SyntaxFunctionCall : public SyntaxElement {

    public:
        SyntaxFunctionCall(const std::string &n, std::list<SyntaxLabeledExpr*>* args);                                  //!< Standard function
        SyntaxFunctionCall(SyntaxVariable*  var, const std::string &n, std::list<SyntaxLabeledExpr*>* args);            //!< Member function
        SyntaxFunctionCall(const SyntaxFunctionCall& x);                                                                //!< Copy constructor
	    
        virtual                            ~SyntaxFunctionCall();                                                       //!< Destructor

        // Assignment operator
        SyntaxFunctionCall&                 operator=(const SyntaxFunctionCall& x);                                     //!< Assignment operator

        // Basic utility functions
        SyntaxFunctionCall*                 clone() const;                                                              //!< Clone object
        
        // Regular functions
        const std::string&                  getFunctionName(void) const { return function_name; }                       //!< Get function name
        RevPtr<RevVariable>                 evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value
        void                                setBaseVariable(SyntaxElement* var) { base_variable = var; }                //!< Set base variable
        std::pair<int,int>                  pipeAddArgPlaceholder(SyntaxElement* piped_arg);                            //!< Add first argument from pipe
        void                                pipeAddArg(SyntaxElement* piped_arg);                                       //!< Add first argument from pipe
        bool                                isConstExpression(void) const;                                              //!< Is this a constant expression?
        bool                                isFunctionSafe(const Environment&       env,
                                                           std::set<std::string>&   localVars) const;                   //!< Is this element safe in a function?
 

    protected:
        
        std::list<SyntaxLabeledExpr*>*      arguments;                                                                  //!< The arguments passed to the function
        std::string                         function_name;                                                               //!< The name of the function
        SyntaxElement*                      base_variable;                                                               //!< Variable holding member function

    };
    
}

#endif

