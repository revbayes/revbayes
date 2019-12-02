/**
 * @file
 * This file contains the declaration of the RevLanguage LG function, which
 * is used to created deterministic variable associated with the fnLg function.
 *
 * @brief Declaration and implementation of RlFunc_lg
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_lg.h 1406 2012-04-20 02:06:14Z awright $
 */

#ifndef Func_lg_H
#define Func_lg_H

#include "Procedure.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_lg :  public Procedure {
        
    public:
        Func_lg( void );
        
        // Basic utility functions
        Func_lg*                                        clone(void) const;                                          //!< Clone the object
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
