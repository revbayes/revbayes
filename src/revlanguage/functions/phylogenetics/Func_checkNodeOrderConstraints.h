/**
 * @file
 * This file contains the declaration of the RevLanguage checkNodeOrderConstraints function, which
 * is used to check whether a given chronogram agrees with some node order relationships.
 *
 * @brief Declaration and implementation of Func_checkNodeOrderConstraints
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2015-11-23 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 */


#ifndef Func_checkNodeOrderConstraints_H
#define Func_checkNodeOrderConstraints_H

#include "RlBoolean.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    class Func_checkNodeOrderConstraints :  public TypedFunction<RlBoolean> {
        
    public:
        Func_checkNodeOrderConstraints( void );
        
        // Basic utility functions
        Func_checkNodeOrderConstraints*                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::Boolean>*               createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
