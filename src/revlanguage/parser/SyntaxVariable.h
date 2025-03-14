#ifndef SyntaxVariable_H
#define SyntaxVariable_H

#include "ModelVector.h"
#include "Natural.h"
#include "SyntaxElement.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    class SyntaxFunctionCall;
    
    
    /**
     * This is the class used to hold variables in the syntax tree.
     *
     * We store the identifier and the base variable
     * here so that we can wrap these things into a DAG node expression
     * if needed.
     *
     * The variable class uses three different functions to evaluate its content.
     * If the variable is part of a left-hand side expression, it is evaluated
     * using evaluateLHSContent(). If it is part of a dynamic rhs expression,
     * it is evaluated using evaluateContent(dynamic=true), and if it is in a static
     * rhs expression, it is evaluated using evaluateContent().
     */
    class SyntaxVariable : public SyntaxElement {
        
    public:
        SyntaxVariable(const std::string &n);                                                                                       //!< Global variable
        SyntaxVariable(const std::string &n, const std::vector<std::string> &ns);                                                                                       //!< Global variable
        
        virtual                            ~SyntaxVariable(void);                                                                   //!< Destructor deletes variable, identifier and index
        
        // Basic utility functions
        SyntaxVariable*                     clone(void) const;                                                                      //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>                 evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);           //!< Get semantic rhs value
        RevPtr<RevVariable>                 evaluateLHSContent(const std::shared_ptr<Environment>& env, const std::string& varType);//!< Get semantic lhs value
        const std::string&                  getIdentifier(void) { return identifier; }                                              //!< Get identifier
        std::string                         getFullName(Environment& env) const;                                                    //!< Get full name, with indices and base obj
        bool                                isFunctionSafe(const Environment&       env,
                                                           std::set<std::string>&   localVars) const;                               //!< Is this element safe in a function?
        bool                                retrievesExternVar(const Environment&       env,
                                                               std::set<std::string>&   localVars,
                                                               bool                     inLHS) const;                               //!< Does this element retrieve an external variable?
        
    protected:

        std::string                         identifier;                                                                             //!< The name of the variable, if identified by name
        std::vector<std::string>            namespaces;
    };
    
}

#endif


