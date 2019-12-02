//
//  Func_DECCladoProbs.h
//  revbayes-proj
//
//  Created by Michael Landis on 1/19/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef Func_DECCladoProbs_H
#define Func_DECCladoProbs_H

#include "RlCladogeneticProbabilityMatrix.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_DECCladoProbs : public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        Func_DECCladoProbs( void );
        
        // Basic utility functions
        Func_DECCladoProbs*                                                clone(void) const;                                      //!< Clone the object
        static const std::string&                                       getClassType(void);                                     //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                 //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                            //!< Get the primary name of the function in Rev
        std::vector<std::string>                                        getFunctionNameAliases( void ) const;
        const TypeSpec&                                                 getTypeSpec(void) const;                                //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix>*         createFunction(void) const;                             //!< Create internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                           //!< Get argument rules
        
    };
    
}


#endif /* defined(__revbayes_proj__Func_DECCladoProbs__) */
