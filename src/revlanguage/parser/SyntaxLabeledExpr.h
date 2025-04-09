#ifndef SyntaxLabeledExpr_H
#define SyntaxLabeledExpr_H

#include "SyntaxElement.h"

#include <iostream>
#include <list>

namespace RevLanguage {

    /**
     * @brief Labeled expression syntax element
     *
     * This syntax element class holds labeled expressions, used
     * in passing arguments to function calls. The member variables
     * are the label and the expression.
     *
     * Labeled expressions occur in function calls like
     *
     *    foo( x = 0.2 )
     *
     * where they correspond to the 'x = 0.2' part. The label may
     * be missing, as it is in the function call
     *
     *    foo( 0.2 )
     *
     * This element just stores the label and expression, and gives
     * them out when the function call syntax element needs them in
     * constructing the function call.s
     */
    class SyntaxLabeledExpr : public SyntaxElement {

    public:
        SyntaxLabeledExpr(const std::string &lbl, SyntaxElement* expr);                                             //!< Constructor
        SyntaxLabeledExpr(const SyntaxLabeledExpr& x);                                                              //!< Copy constructor
	    
        virtual                    ~SyntaxLabeledExpr();                                                            //!< Destructor

        // Assignment operator
        SyntaxLabeledExpr&          operator=(const SyntaxLabeledExpr& x);                                          //!< Assignment operator

        // Basic utility functions
        SyntaxLabeledExpr*          clone(void) const;                                                              //!< Clone object
        
        // Regular functions
        const SyntaxElement&        getExpression(void) const { return *expression; }                               //!< Return expression
        SyntaxElement&              getExpression(void) { return *expression; }                                     //!< Return expression
        const std::string&          getLabel() const { return label; }                                              //!< Return label
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value
        bool                        isConstExpression(void) const;                                                  //!< Is the expression constant?
        bool                        isFunctionSafe(const Environment&       env,
                                                   std::set<std::string>&   localVars) const;                       //!< Is this element safe in a function?
        bool                        retrievesExternVar(const Environment&       env,
                                                       std::set<std::string>&   localVars,
                                                       bool                     inLHS) const;                       //!< Does this element retrieve an external variable?

    protected:
        
        std::string                 label;                                                                          //!< The label of the argument
        SyntaxElement*              expression;                                                                     //!< The expression for the argument value
    
    };
    
}

#endif

