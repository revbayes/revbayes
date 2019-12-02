//
//  RlCharacterHistoryRateModifier.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "RlCharacterHistoryRateModifier.h"
#include "CharacterHistoryRateModifier.h"
#include "ArgumentRule.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RlBoolean.h"

using namespace RevLanguage;

CharacterHistoryRateModifier::CharacterHistoryRateModifier(void) : ModelObject<RevBayesCore::CharacterHistoryRateModifier>() {
    
}


CharacterHistoryRateModifier::CharacterHistoryRateModifier( const RevBayesCore::CharacterHistoryRateModifier &v) : ModelObject<RevBayesCore::CharacterHistoryRateModifier>( v.clone() )
{
    
}

CharacterHistoryRateModifier::CharacterHistoryRateModifier( RevBayesCore::CharacterHistoryRateModifier *v) : ModelObject<RevBayesCore::CharacterHistoryRateModifier>( v ) {
    
}


CharacterHistoryRateModifier::CharacterHistoryRateModifier( RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryRateModifier> *m) : ModelObject<RevBayesCore::CharacterHistoryRateModifier>( m ) {
    
}


CharacterHistoryRateModifier* CharacterHistoryRateModifier::clone() const {
    return new CharacterHistoryRateModifier( *this );
}


/* Get Rev type of object */
const std::string& CharacterHistoryRateModifier::getClassType(void) {
    
    static std::string rev_type = "CharacterHistoryRateModifier";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& CharacterHistoryRateModifier::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& CharacterHistoryRateModifier::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
