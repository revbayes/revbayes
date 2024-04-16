//
//  Func_vectorFlatten.h
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__Func_vectorFlatten__
#define __revbayes_proj__Func_vectorFlatten__

#include "ModelVector.h"
#include "Real.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    class Func_vectorFlatten :  public TypedFunction<ModelVector<Real> > {
        
    public:
        Func_vectorFlatten();
        
        // Basic utility functions
        Func_vectorFlatten*                                                     clone(void) const;                                          //!< Clone the object
        static const std::string&                                               getClassType(void);                                         //!< Get Rev type (static)
        static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get Rev type spec (static)
        std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get Rev type spec (from instance)
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector< double> >*         createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}


#endif /* defined(__revbayes_proj__Func_vectorFlatten__) */
