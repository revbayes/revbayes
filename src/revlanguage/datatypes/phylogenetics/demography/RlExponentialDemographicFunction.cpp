#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ExponentialDemographicFunction.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlString.h"
#include "RlExponentialDemographicFunction.h"
#include "TypeSpec.h"


using namespace RevLanguage;

ExponentialDemographicFunction::ExponentialDemographicFunction(void) : DemographicFunction()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
ExponentialDemographicFunction* ExponentialDemographicFunction::clone(void) const
{
    
    return new ExponentialDemographicFunction(*this);
}


void ExponentialDemographicFunction::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate space for a new ExponentialDemographicFunction object
    RevBayesCore::TypedDagNode<double> *t0 = static_cast<const RealPos &>( time_recent->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double> *t1 = static_cast<const RealPos &>( time_ancient->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double> *N0 = static_cast<const RealPos &>( theta_recent->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double> *N1 = static_cast<const RealPos &>( theta_ancient->getRevObject() ).getDagNode();

    value = new RevBayesCore::ExponentialDemographicFunction( N0, N1, t0, t1 );
}


/** Get Rev type of object */
const std::string& ExponentialDemographicFunction::getClassType(void)
{
    
    static std::string rev_type = "ExponentialDemographicFunction";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& ExponentialDemographicFunction::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( DemographicFunction::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string ExponentialDemographicFunction::getDemographicFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Exponential";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& ExponentialDemographicFunction::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        member_rules.push_back( new ArgumentRule("N0", RealPos::getClassTypeSpec(), "The population size at the beginning of the growth period (towards the present).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("N1", RealPos::getClassTypeSpec(), "The population size at the end of the growth period (towards the past).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("t0", RealPos::getClassTypeSpec(), "The time of the beginning of the growth period (towards the present).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("t1", RealPos::getClassTypeSpec(), "The time of the end of the growth period (towards the past).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return member_rules;
}

/** Get type spec */
const TypeSpec& ExponentialDemographicFunction::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void ExponentialDemographicFunction::printValue(std::ostream &o) const
{
    
    o << "ExponentialDemographicFunction";
}


/** Set a member variable */
void ExponentialDemographicFunction::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "N0" )
    {
        theta_recent = var;
    }
    else if ( name == "N1" )
    {
        theta_ancient = var;
    }
    else if ( name == "t0" )
    {
        time_recent = var;
    }
    else if ( name == "t1" )
    {
        time_ancient = var;
    }
    else
    {
        DemographicFunction::setConstParameter(name, var);
    }
}
