/**
 * @file
 * This file contains the declaration of the RevLanguage gtr function, which
 * is used to created deterministic variable associated with the gtr function.
 *
 * @brief Declaration and implementation of Func_gtr
 *
 * (c) Copyright 2009- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */


#ifndef Func_dppMeanFromConc_H
#define Func_dppMeanFromConc_H

#include "Procedure.h"

#include <string>

namespace RevLanguage {
    
    class Func_dppMeanFromConc :  public Procedure {
        
    public:
        Func_dppMeanFromConc( void );
        
        // Basic utility functions
        Func_dppMeanFromConc*                           clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance

        // Function functions you have to override
        RevPtr<RevVariable>                             execute(void);                                              //!< Execute function
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                                 getReturnType(void) const;                                  //!< Get type of return value
        
    };
    
}

#endif
