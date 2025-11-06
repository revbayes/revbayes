#ifndef SyntaxForLoop_H
#define SyntaxForLoop_H

#include "SyntaxElement.h"

#include <iostream>
#include <list>

namespace RevLanguage {

    /**
     * @brief For loop syntax element
     *
     * This syntax element is used to hold for loops, more specifically
     * the definition of the loop, that is the part
     *
     *    for ( i in 1:10 )
     *
     * The entire loop construct is held by the statement syntax element.
     *
     * Here, we store the loop variable name and the in-expression. We also
     * help the statement syntax element execute the for loop by providing
     * functions that initialize the loop, finalize the loop, get the next
     * loop state, and test whether the loop has finished.
     *
     * When the loop is initialized, we put the loop variable in the
     * execution environment as a control variable, to make sure that it
     * is not included in DAGs if it appears in dynamic expressions. We
     * then assign values to it from the variable resulting from execution
     * of the in-expression. This variable needs to be a container with at
     * least one dimension (a vector), and we get consecutive slices out of
     * it from its first dimension. For a vector, this simply means that the
     * loop variable takes on each of the values of the vector in turn.
     *
     * Like in R, loops do not open up new local environment. All statements
     * in the for loop are executed in the outer environment, and the loop
     * variable remains there after the loop finishes, as in R.
     */
    class SyntaxForLoop : public SyntaxElement {

    public:
        SyntaxForLoop(const std::string &identifier, SyntaxElement* inExpr);                       //!< Standard constructor
        SyntaxForLoop(const SyntaxForLoop& x);                                                     //!< Copy constructor
	    
        virtual                     ~SyntaxForLoop();                                               //!< Destructor

        // Assignment operator
        SyntaxForLoop&              operator=(const SyntaxForLoop& x);                              //!< Assignment operator
        
        // Basic utility functions
        SyntaxElement*              clone() const;                                                  //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value
        void                        finalizeLoop(void);                                             //!< Finalize loop
        const std::string&          getIndexVarName(void) const;                                    //!< Get the name of the index variable
        void                        getNextLoopState(void);                                         //!< Get next state of loop
        bool                        isFinished() const;                                             //!< Have we iterated over the whole loop?
        void                        initializeLoop(const std::shared_ptr<Environment>& env);        //!< Initialize loop
        bool                        isFunctionSafe(const Environment&       env,
                                                   std::set<std::string>&   localVars) const;       //!< Is this element safe in a function?

    protected:
        
        std::string                 varName;                                                        //!< The name of the loop variable
        SyntaxElement*              inExpression;                                                   //!< The in expression (a vector of values)
        Container*                  stateSpace;                                                     //!< Vector result of 'in' expression
        size_t                      nextIndex;                                                      //!< Next element in vector
        RevPtr<RevVariable>         loopVariable;                                                   //!< Smart pointer to the loop variable in the environment

    };
    
}

#endif

