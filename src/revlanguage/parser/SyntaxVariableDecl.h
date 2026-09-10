#ifndef SyntaxVariableDecl_H
#define SyntaxVariableDecl_H

#include "SyntaxElement.h"
#include "SyntaxVariable.h"

#include <iostream>
#include <list>
#include <string>

namespace RevLanguage {

    /**
     * @brief Variable declarations
     *
     * This syntax element class deals with variable declarations.
     * They are used to define the type of variables in the workspace,
     * and the lengths (size) and dimensions of generic containers. For
     * instance,
     * 
     * Real[8] x
     *
     * would declare and create a generic vector of Real numbers of
     * length 8, containing default Real elements.
     *
     * Similar expressions are used to define types of variables in
     * function definitions and class definitions. Such expressions are
     * handled by SyntaxFormal. The principal difference between a
     * formal and a variable declaration is that the former do no
     * support lengths specifications, whereas variable declarations
     * do.
     *
     * Note that it is possible to leave out the the lengths specification.
     * For instance,
     *
     * Real[] x
     *
     * would specify a generic vector of Real numbers, which will be
     * initialized as an empty vector. This is perhaps going to be the
     * typical usage of variable declarations.
     */
    class SyntaxVariableDecl : public SyntaxElement {
    
    public:
        SyntaxVariableDecl(const std::string &typeName, std::list<SyntaxElement*>* lengths, const std::string &varName);    //!< Basic constructor
        SyntaxVariableDecl(const SyntaxVariableDecl& x);                                                                    //!< Copy constructor
        
        virtual                    ~SyntaxVariableDecl();                                                                   //!< Virtual destructor
        
        // Assignment operator
        SyntaxVariableDecl&         operator=(const SyntaxVariableDecl& x);                     //!< Assignment operator
        
        // Basic utility functions
        SyntaxVariableDecl*         clone() const;                                              //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);      //!< Get semantic value
        bool                        isFunctionSafe(const Environment&       env,
                                                   std::set<std::string>&   localVars) const;   //!< Is this element safe in a function?
        
    protected:
        
        std::string                 elementTypeName;                                            //!< Element type of the variable
        std::list<SyntaxElement*>*  lengthExpr;                                                 //!< Lengths in different dimensions
        std::string                 variableName;                                               //!< Variable name
    
    };
    
}

#endif

