//
//  Func_shortestDistance.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/15/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef Func_shortestDistance_h
#define Func_shortestDistance_h

#include "ModelVector.h"
#include "RbVector.h"
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
    class Func_shortestDistance :  public TypedFunction< ModelVector<ModelVector<RealPos> > > {
        
    public:
        Func_shortestDistance( void );
        
        // Basic utility functions
        Func_shortestDistance*                              clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}




#endif /* Func_shortestDistance_h */
