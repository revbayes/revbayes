//
//  RlCharacterHistoryRateModifier.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef RlCharacterHistoryRateModifier_h
#define RlCharacterHistoryRateModifier_h

#include "ModelObject.h"
#include "CharacterHistoryRateModifier.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class CharacterHistoryRateModifier : public ModelObject<RevBayesCore::CharacterHistoryRateModifier> {
        
    public:
        
        CharacterHistoryRateModifier(void);                                                                                //!< Default constructor
        CharacterHistoryRateModifier(const RevBayesCore::CharacterHistoryRateModifier& m);                                 //!< Default constructor
        CharacterHistoryRateModifier(RevBayesCore::CharacterHistoryRateModifier *m);                                       //!< Default constructor
        CharacterHistoryRateModifier(RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryRateModifier> *d);           //!< Default constructor
        
        // Basic utility functions
        CharacterHistoryRateModifier*       clone(void) const;                                                          //!< Clone object
        static const std::string&           getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                     //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                    //!< Get language type of the object
        
    };
    
}

#endif /* RlCharacterHistoryRateModifier_h */
