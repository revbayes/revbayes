//
//  Func_EpochCladoProbs.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/24/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Func_EpochCladoProbs_h
#define Func_EpochCladoProbs_h

#include "RlCladogeneticProbabilityMatrix.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_EpochCladoProbs : public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        Func_EpochCladoProbs( void );
        
        // Basic utility functions
        Func_EpochCladoProbs*                                           clone(void) const;                                      //!< Clone the object
        static const std::string&                                       getClassType(void);                                     //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                 //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                            //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                                //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix>*         createFunction(void) const;                             //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                           //!< Get argument rules
        
    };
    
}

#endif /* Func_EpochCladoProbs_h */
