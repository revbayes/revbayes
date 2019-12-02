#include "Probability.h"
#include "RealPos.h"
#include "RbException.h"
#include "RlString.h"
#include "RbUtil.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

#include <sstream>

using namespace RevLanguage;

/** Default constructor */
Probability::Probability( void ) : RealPos( 1.0 )
{

}


/** Construct from double */
Probability::Probability( double x ) : RealPos( x )
{
    
    if ( x < 0.0 || x > 1.0)
    {
        throw RbException( "Creation of " + getClassType() + " with value x=" + x + " outside standard probabilities (0,1)");
    }
    
}


/** Construct from double */
Probability::Probability( RevBayesCore::TypedDagNode<double> *x ) : RealPos( x )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Probability* Probability::clone( void ) const
{
    
	return new Probability( *this );
}

RevObject* Probability::convertTo( const TypeSpec& type ) const
{
    
    if ( type == RealPos::getClassTypeSpec() )
    {
        return new RealPos(dag_node->getValue());
    }
    
    return RealPos::convertTo( type );
}


/** Get Rev type of object */
const std::string& Probability::getClassType(void)
{
    
    static std::string rev_type = "Probability";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Probability::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RealPos::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& Probability::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Is convertible to type? */
double Probability::isConvertibleTo( const TypeSpec& type, bool once ) const
{
    
    if ( type == RealPos::getClassTypeSpec() )
    {
        return 0.1;
    }
    else
    {
        double tmp = RealPos::isConvertibleTo(type, once);
        return ( (tmp == -1.0) ? -1.0 : (tmp+0.1));
    }
}

