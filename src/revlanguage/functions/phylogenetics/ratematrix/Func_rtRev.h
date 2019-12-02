/**
 * @file
 * This file contains the declaration of the RevLanguage RtRev function, which
 * is used to created deterministic variable associated with the rtRev function.
 *
 * @brief Declaration and implementation of Func_rtRev
 *
 * (c) Copyright 2009- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */


#ifndef Func_rtRev_H
#define Func_rtRev_H

#include "Procedure.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_rtRev :  public Procedure {
        
    public:
        Func_rtRev( void );
        
        // Basic utility functions
        Func_rtRev*                                     clone(void) const;                                          //!< Clone the object
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
