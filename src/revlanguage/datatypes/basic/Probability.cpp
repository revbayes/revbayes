#include "Probability.h"

#include <string>

#include "RealPos.h"
#include "RbException.h"
#include "ReferenceFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevLanguage { class RevObject; }

using namespace RevLanguage;

/** Default constructor */
Probability::Probability( void ) : RealPos( 1.0 )
{

}


/** Construct from double */
Probability::Probability( double x ) : RealPos( x )
{
    validate();
}


/** Construct from double */
Probability::Probability( RevBayesCore::TypedDagNode<double> *x ) : RealPos( x )
{
    validate();
}

void Probability::validate(double x) const
{
    if (x < 0 or x > 1)
	throw RbException() << std::setprecision(17) << "Creation of " << getClassType() << " with value x=" << x << " outside standard probabilities [0,1]";
}

void Probability::validate() const
{
    return validate(dag_node->getValue());
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


void Probability::constructInternalObject( void )
{
  // we free the memory first
    if ( dag_node != NULL )
    {
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
    }

    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( real->getRevObject() ).getDagNode();
    RevBayesCore::ReferenceFunction<double>* f = new RevBayesCore::ReferenceFunction<double>( x );


    dag_node = new RevBayesCore::DeterministicNode<double>("", f);
    dag_node->incrementReferenceCount();

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


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Probability::getConstructorFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "Probability";

    return f_name;
}


/** Return member rules (no members) */
const MemberRules& Probability::getParameterRules(void) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
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


/** Set a member variable */
void Probability::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "x" )
    {
        real = var;
	RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( real->getRevObject() ).getDagNode();
	validate(x->getValue());
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
