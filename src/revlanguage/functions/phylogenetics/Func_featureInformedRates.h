//
//  Func_featureInformedRates.hpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/12/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#ifndef Func_featureInformedRates_hpp
#define Func_featureInformedRates_hpp


#include <string>
#include <iosfwd>
#include <vector>

#include "RlCladogeneticProbabilityMatrix.h"
#include "RlTypedFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_featureInformedRates : public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_featureInformedRates( void );
        
        // Basic utility functions
        Func_featureInformedRates*                                                clone(void) const;                                      //!< Clone the object
        static const std::string&                                       getClassType(void);                                     //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                 //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                            //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >*         createFunction(void) const;                             //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                           //!< Get argument rules
        
    };
    
}
#endif /* Func_featureInformedRates_hpp */
