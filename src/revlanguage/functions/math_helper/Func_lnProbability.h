//
//  Func_lnProbability.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/12/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__Func_lnProbability__
#define __revbayes_proj__Func_lnProbability__


#include "Probability.h"
#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage conversion function between a Real and a Probability.
     *
     * This Rev functions converts between Rev types from Real to Probability.
     * The Rev conversion functions only checks during instantiation if the conversion
     * is possible. We then simply return the same DAG node but in a different Rev wrapper.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-31, version 1.0
     *
     */
    class Func_lnProbability : public TypedFunction<Real> {
        
    public:
        Func_lnProbability( void );
        
        // Basic utility functions
        Func_lnProbability*                             clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif /* defined(__revbayes_proj__Func_lnProbability__) */
