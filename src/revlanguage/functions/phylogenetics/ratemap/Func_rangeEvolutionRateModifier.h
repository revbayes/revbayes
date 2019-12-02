//
//  Func_rangeEvolutionRateModifier.hpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef Func_rangeEvolutionRateModifier_hpp
#define Func_rangeEvolutionRateModifier_hpp

#include "RlCharacterHistoryRateModifier.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_rangeEvolutionRateModifier : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        Func_rangeEvolutionRateModifier( void );
        
        // Basic utility functions
        Func_rangeEvolutionRateModifier*                                                  clone(void) const;                       //!< Clone the object
        static const std::string&                                                   getClassType(void);                      //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                  //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;             //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                 //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >*  createFunction(void) const;              //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;            //!< Get argument rules
        
    };
    
}


#endif /* Func_rangeEvolutionRateModifier_hpp */
