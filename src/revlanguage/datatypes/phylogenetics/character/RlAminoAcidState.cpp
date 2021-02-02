#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "RlAminoAcidState.h"
#include "TypeSpec.h"
#include "AminoAcidState.h"
#include "CharacterState.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
AminoAcidState::AminoAcidState(void) : ModelObject<RevBayesCore::AminoAcidState>()
{
    
}

/** Construct from bool */
AminoAcidState::AminoAcidState(const RevBayesCore::AminoAcidState &d) : ModelObject<RevBayesCore::AminoAcidState>( new RevBayesCore::AminoAcidState(d) ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
AminoAcidState* AminoAcidState::clone(void) const
{
    
	return new AminoAcidState(*this);
}


/** Get Rev type of object */
const std::string& AminoAcidState::getClassType(void) { 
    
    static std::string rev_type = "AminoAcid";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& AminoAcidState::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& AminoAcidState::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

