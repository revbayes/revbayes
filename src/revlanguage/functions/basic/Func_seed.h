/**
 * @file
 * This file contains the declaration and implementation
 * of the Func_seed, which is used to set the seed of the random number generator.
 *
 * @brief Declaration and implementation of Func_seed
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-12 10:25:58 +0200 (Tue, 12 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: Func_seed.h 1626 2012-06-12 08:25:58Z hoehna $
 */

#ifndef Func_seed_H
#define Func_seed_H

#include "Procedure.h"

#include <string>

namespace RevLanguage {
    
    class Func_seed : public Procedure {
        
    public:
        Func_seed();
        
        // Basic utility functions
        Func_seed*                  clone(void) const;                                          //!< Clone the object
        static const std::string&   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&             getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        const ArgumentRules&        getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&             getReturnType(void) const;                                  //!< Get type of return value
        
        
        RevPtr<RevVariable>         execute(void);                                              //!< Execute function
        
    };
    
}


#endif


