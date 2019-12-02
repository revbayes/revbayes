/**
 * @file
 * This file contains the declaration of ConstructorFunction, which is used
 * for functions that construct objects using the RevLanguage.
 *
 * @brief Declaration of ConstructorFunction
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-01 14:40:37 +0200 (Fri, 01 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: ConstructorFunction.h 1602 2012-06-01 12:40:37Z hoehna $
 */

#ifndef ConstructorFunction_H
#define ConstructorFunction_H

#include "Procedure.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace RevLanguage {

class ConstructorFunction : public Procedure {
    
public:
    ConstructorFunction(RevObject *obj);                                                                    //!< Object constructor
    ConstructorFunction(const ConstructorFunction& obj);                                                    //!< Copy constructor
    virtual ~ConstructorFunction(void);
    
    // overloaded operators
    ConstructorFunction&                    operator=(const ConstructorFunction& c);
    
    // Basic utility functions
    ConstructorFunction*                    clone(void) const;                                              //!< Clone the object
    static const std::string&               getClassType(void);                                             //!< Get Rev type
    static const TypeSpec&                  getClassTypeSpec(void);                                         //!< Get class type spec
    std::string                             getFunctionName(void) const;                                    //!< Get the primary name of the function in Rev
    std::vector<std::string>                getFunctionNameAliases(void) const;                             //!< Get the aliases of the name of the function in Rev
    const TypeSpec&                         getTypeSpec(void) const;                                        //!< Get language type of the object
    
    // Regular functions
    RevPtr<RevVariable>                     execute(void);                                                  //!< Execute the function. This is the function one has to overwrite for single return values.
    const ArgumentRules&                    getArgumentRules(void) const;                                   //!< Get argument rules
    const TypeSpec&                         getReturnType(void) const;                                      //!< Get type of return value
    RevObject*                              getRevObject(void) { return templateObject; }                   //!< Get type of return value
    
protected:
    
    const ArgumentRules*                    argRules;                                                       //!< Member rules converted to reference rules
    RevObject*                              templateObject;                                                 //!< The template object
    
};
    
}

#endif
