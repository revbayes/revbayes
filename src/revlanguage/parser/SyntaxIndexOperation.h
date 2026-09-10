#ifndef SyntaxIndexOperation_H
#define SyntaxIndexOperation_H

#include "SyntaxElement.h"

#include <iostream>
#include <list>

namespace RevLanguage {

    class SyntaxFunctionCall;
    class RevVariableSlot;

    /**
     * This is the class used to hold variables in the syntax tree.
     *
     * We store the identifier, the index vector and the base variable
     * here so that we can wrap these things into a DAG node expression
     * if needed.
     *
     */

    class SyntaxIndexOperation : public SyntaxElement {

    public:
        SyntaxIndexOperation(SyntaxElement* var, SyntaxElement* indx);                                                              //!< Standard constructor
        SyntaxIndexOperation(const SyntaxIndexOperation& x);                                                                        //!< Copy constructor

        virtual                            ~SyntaxIndexOperation(void);                                                             //!< Destructor deletes variable, identifier and index

        // Assignment operator
        SyntaxIndexOperation&               operator=(const SyntaxIndexOperation& x);                                               //!< Assignment operator

        // Basic utility functions
        SyntaxIndexOperation*               clone(void) const;                                                                      //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>                 evaluateLHSContent(const std::shared_ptr<Environment>& env, const std::string& varType);   //!< Get semantic lhs value
        RevPtr<RevVariable>                 evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);           //!< Get semantic value
        SyntaxElement*                      getBaseVariable(void);                                                                  //!< Get the base variable for this expression
        void                                updateVariable(Environment& env, const std::string &n);                                 //!< Update the composite variables
        
    protected:
        SyntaxElement*                      index;                                                                                  //!< Vector of int indices to variable element
        SyntaxElement*                      base_variable;                                                                           //!< Base variable (pointing to a composite node)
    };

}

#endif

