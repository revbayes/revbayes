//
//  Func_empiricalQuantile.hpp
//  revbayes
//
//  Created by Lyndon Coghill on 9/22/16.
//  Copyright © 2016 Lyndon Coghill. All rights reserved.
//

#ifndef Func_empiricalQuantile_H
#define Func_empiricalQuantile_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the arithmetic empiricalQuantile function.
     *
     * The RevLanguage wrapper of the empiricalQuantile function connects
     * the variables/parameters of the function and creates the internal empiricalQuantileFunction object.
     * Please read the empiricalQuantileFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_empiricalQuantile :  public TypedFunction<Real> {
        
    public:
        Func_empiricalQuantile( void );
        
        // Basic utility functions
        Func_empiricalQuantile*                                      clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
