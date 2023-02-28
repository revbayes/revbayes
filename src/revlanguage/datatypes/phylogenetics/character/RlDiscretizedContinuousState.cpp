#include "RlDiscretizedContinuousState.h"

#include "ConstantNode.h"
#include "TypeSpec.h"
#include "RevObject.h"

using namespace RevLanguage;

/** Default constructor */
DiscretizedContinuousState::DiscretizedContinuousState(void) : ModelObject<RevBayesCore::DiscretizedContinuousState>()
{
    
}

/** Construct from bool */
DiscretizedContinuousState::DiscretizedContinuousState(const RevBayesCore::DiscretizedContinuousState &d) : ModelObject<RevBayesCore::DiscretizedContinuousState>( new RevBayesCore::DiscretizedContinuousState(d) )
{
    
}

DiscretizedContinuousState::DiscretizedContinuousState( RevBayesCore::TypedDagNode<RevBayesCore::DiscretizedContinuousState> *v ) : ModelObject<RevBayesCore::DiscretizedContinuousState>( v )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
DiscretizedContinuousState* DiscretizedContinuousState::clone(void) const {
    
	return new DiscretizedContinuousState(*this);
}



/** Get Rev type of object */
const std::string& DiscretizedContinuousState::getClassType(void) {
    
    static std::string rev_type = "DiscretizedContinuousState";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& DiscretizedContinuousState::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}



/** Get type spec */
const TypeSpec& DiscretizedContinuousState::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}




