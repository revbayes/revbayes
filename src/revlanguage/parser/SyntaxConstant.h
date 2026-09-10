#ifndef SyntaxConstant_H
#define SyntaxConstant_H

#include "SyntaxElement.h"

#include <iostream>


namespace RevLanguage {

    /**
     * This is the class used to hold constants in the syntax tree.
     *
     * The result is never going to change so we can set the value
     * already in the constructor.
     *
     */
    class SyntaxConstant : public SyntaxElement {

    public:
                                                SyntaxConstant(RevObject* val);                         //!< Constructor from value
                                                SyntaxConstant(const SyntaxConstant& sc);               //!< Copy constructor
	    virtual                                ~SyntaxConstant();                                       //!< Destructor deletes value

        // Assignment operator
        SyntaxConstant&                         operator=(const SyntaxConstant& x);                     //!< Assignment operator

        // Basic utility functions
        SyntaxConstant*                         clone(void) const;                                      //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>                     evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);  //!< Get semantic value
        bool                                    isConstExpression(void) const;                          //!< Is the expression constant?

    protected:
        
        RevObject*                              value;                                                  //!< The constant value
    
    };
    
}

#endif

