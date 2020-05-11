//
//  Func_chromosomesCladoProbs.h
//
//  Created by will freyman on 12/14/15.
//  Copyright (c) 2015 will freyman. All rights reserved.
//

#ifndef Func_chromosomesCladoProbs_H
#define Func_chromosomesCladoProbs_H


#include <string>
#include <iosfwd>
#include <vector>

#include "RlCladogeneticProbabilityMatrix.h"
#include "RlTypedFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_chromosomesCladoProbs : public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        Func_chromosomesCladoProbs( void );
        
        // Basic utility functions
        Func_chromosomesCladoProbs*                                                     clone(void) const;                                      //!< Clone the object
        static const std::string&                                                       getClassType(void);                                     //!< Get Rev type
        static const TypeSpec&                                                          getClassTypeSpec(void);                                 //!< Get class type spec
        std::string                                                                     getFunctionName(void) const;                            //!< Get the primary name of the function in Rev
        const TypeSpec&                                                                 getTypeSpec(void) const;                                //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix>*      createFunction(void) const;                             //!< Create internal function object
        const ArgumentRules&                                                            getArgumentRules(void) const;                           //!< Get argument rules
        
    };
    
}




#endif
