//
//  Func_hostSwitchRateModifier.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/24/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef Func_hostSwitchRateModifier_h
#define Func_hostSwitchRateModifier_h

#include <string>
#include <iosfwd>
#include <vector>

#include "RlCharacterHistoryRateModifier.h"
#include "RlTypedFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevLanguage {
class ArgumentRules;
class TypeSpec;
    
    class Func_hostSwitchRateModifier : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        Func_hostSwitchRateModifier( void );
        
        // Basic utility functions
        Func_hostSwitchRateModifier*                                      clone(void) const;                       //!< Clone the object
        static const std::string&                                                   getClassType(void);                      //!< Get Rev type
        static const TypeSpec&                                                      getClassTypeSpec(void);                  //!< Get class type spec
        std::string                                                                 getFunctionName(void) const;             //!< Get the primary name of the function in Rev
        const TypeSpec&                                                             getTypeSpec(void) const;                 //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >*  createFunction(void) const;              //!< Create internal function object
        const ArgumentRules&                                                        getArgumentRules(void) const;            //!< Get argument rules
        
    };
    
}
#endif /* Func_hostSwitchRateModifier_h */
