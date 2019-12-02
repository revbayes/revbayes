#include "RlNaturalNumbersState.h"
#include "ConstantNode.h"
#include "RbUtil.h"
#include "TypeSpec.h"
#include <sstream>

using namespace RevLanguage;

/** Default constructor */
NaturalNumbersState::NaturalNumbersState(void) : ModelObject<RevBayesCore::NaturalNumbersState>()
{
    
}

/** Construct from bool */
NaturalNumbersState::NaturalNumbersState(const RevBayesCore::NaturalNumbersState &d) : ModelObject<RevBayesCore::NaturalNumbersState>( new RevBayesCore::NaturalNumbersState(d) )
{
    
}

NaturalNumbersState::NaturalNumbersState( RevBayesCore::TypedDagNode<RevBayesCore::NaturalNumbersState> *v ) : ModelObject<RevBayesCore::NaturalNumbersState>( v )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
NaturalNumbersState* NaturalNumbersState::clone(void) const {
    
	return new NaturalNumbersState(*this);
}



/** Get Rev type of object */
const std::string& NaturalNumbersState::getClassType(void) {
    
    static std::string rev_type = "NaturalNumbers";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& NaturalNumbersState::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}



/** Get type spec */
const TypeSpec& NaturalNumbersState::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}




