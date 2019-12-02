//
//  Func_geographicalDistance.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Func_geographicalDistance_h
#define Func_geographicalDistance_h

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>


namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the hyperbolic Sine function (sinh()).
     *
     * The RevLanguage wrapper of the hyperbolic function function connects
     * the variables/parameters of the function and creates the internal HyperbolicSineFunction object.
     * Please read the HyperbolicSineFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since 2016-09-26, version 1.0
     *
     */
    class Func_geographicalDistance :  public TypedFunction<RealPos> {
        
    public:
        Func_geographicalDistance( void );
        
        // Basic utility functions
        Func_geographicalDistance*                         clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* Func_geographicalDistance_h */
