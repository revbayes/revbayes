#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "Natural.h"
#include "NaturalNumbersState.h"
#include "RlDiscreteCharacterState.h"
#include "StandardState.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "CharacterState.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
DiscreteCharacterState::DiscreteCharacterState(void) : ModelObject<RevBayesCore::DiscreteCharacterState>()
{
    
}

/** Construct from bool */
DiscreteCharacterState::DiscreteCharacterState(const RevBayesCore::DiscreteCharacterState &d) : ModelObject<RevBayesCore::DiscreteCharacterState>( d.clone() )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
DiscreteCharacterState* DiscreteCharacterState::clone(void) const
{
    
    return new DiscreteCharacterState(*this);
}


/** Convert to type. The caller manages the returned object. */
RevObject* DiscreteCharacterState::convertTo( const TypeSpec& type ) const
{
    
    
    if ( type == Natural::getClassTypeSpec() )
    {
        const RevBayesCore::NaturalNumbersState* ns =dynamic_cast<const RevBayesCore::NaturalNumbersState*>( &getValue() );
        if ( ns != NULL )
        {
            return new Natural( StringUtilities::asIntegerNumber( ns->getStringValue() ) );
        }
        
        const RevBayesCore::StandardState* ss =dynamic_cast<const RevBayesCore::StandardState*>( &getValue() );
        if ( ss != NULL && StringUtilities::isIntegerNumber( ss->getStringValue() ) )
        {
            return new Natural( StringUtilities::asIntegerNumber( ss->getStringValue() ) );
        }
        
    }
    
    return ModelObject<RevBayesCore::DiscreteCharacterState>::convertTo( type );
}



/** Get Rev type of object */
const std::string& DiscreteCharacterState::getClassType(void)
{
    
    static std::string rev_type = "DiscreteCharacter";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& DiscreteCharacterState::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& DiscreteCharacterState::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Is convertible to type? */
double DiscreteCharacterState::isConvertibleTo( const TypeSpec& type, bool once ) const
{
    
    if ( once == true && type == Natural::getClassTypeSpec() )
    {
        const RevBayesCore::NaturalNumbersState* ns =dynamic_cast<const RevBayesCore::NaturalNumbersState*>( &getValue() );
        if ( ns != NULL )
        {
            return 0.1;
        }
        
        const RevBayesCore::StandardState* ss =dynamic_cast<const RevBayesCore::StandardState*>( &getValue() );
        if ( ss != NULL && StringUtilities::isIntegerNumber( ss->getStringValue() ) )
        {
            return 0.1;
        }
        
    }
    
    return ModelObject<RevBayesCore::DiscreteCharacterState>::isConvertibleTo( type, once );
}

