//
//  Func_gamma.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/21/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Func_gamma_h
#define Func_gamma_h

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the gamma function.
     *
     * The RevLanguage wrapper of the gamma function connects
     * the variables/parameters of the function and creates the internal GammaFunction object.
     * Please read the GammaFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2016-10-21, version 1.0
     *
     */
    class Func_gamma : public TypedFunction<RealPos> {
        
    public:
        Func_gamma( void );
        
        // Basic utility functions
        Func_gamma*                                     clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* Func_gamma_h */
