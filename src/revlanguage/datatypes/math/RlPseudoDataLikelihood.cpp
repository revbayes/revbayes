#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "RlPseudoDataLikelihood.h"

using namespace RevLanguage;

/* Default constructor */
PseudoDataLikelihood::PseudoDataLikelihood(void) : ModelObject<RevBayesCore::PseudoDataLikelihood>( new RevBayesCore::PseudoDataLikelihood() )
{
}

PseudoDataLikelihood::PseudoDataLikelihood(const RevBayesCore::PseudoDataLikelihood& from) : ModelObject<RevBayesCore::PseudoDataLikelihood>( new RevBayesCore::PseudoDataLikelihood(from) )
{
}

PseudoDataLikelihood::PseudoDataLikelihood(RevBayesCore::PseudoDataLikelihood* l) : ModelObject<RevBayesCore::PseudoDataLikelihood>( l )
{
}

PseudoDataLikelihood::PseudoDataLikelihood( RevBayesCore::TypedDagNode<RevBayesCore::PseudoDataLikelihood> * l ) : ModelObject<RevBayesCore::PseudoDataLikelihood>( l )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
PseudoDataLikelihood* PseudoDataLikelihood::clone(void) const
{
	return new PseudoDataLikelihood(*this);
}

/** Get Rev type of object */
const std::string& PseudoDataLikelihood::getClassType(void)
{
    
    static std::string rev_type = "PseudoDataLikelihood";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& PseudoDataLikelihood::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get type spec */
const TypeSpec& PseudoDataLikelihood::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


