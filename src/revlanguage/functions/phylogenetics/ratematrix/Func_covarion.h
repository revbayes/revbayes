//
//  Func_covarion.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/4/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef Func_covarion_h
#define Func_covarion_h

#include "RateGenerator.h"
#include "ModelVector.h"
#include "RlRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_covarion : public TypedFunction<RateMatrix> {
        
    public:
        Func_covarion( void );
        
        // Basic utility functions
        Func_covarion*                                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >*     createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* Func_covarion_h */
