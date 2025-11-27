#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_geom.h"
#include "GeometricDistribution.h"
#include "Natural.h"
#include "Probability.h"
#include "StochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_geom::Dist_geom() : TypedDistribution<Natural>() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
Dist_geom* Dist_geom::clone( void ) const 
{
    
    return new Dist_geom(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be 
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the 
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::GeometricDistribution* Dist_geom::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* prob    = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::GeometricDistribution* d      = new RevBayesCore::GeometricDistribution( prob );
    
    return d;
}



/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_geom::getClassType(void) 
{ 
    
    static std::string rev_type = "Dist_geom";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_geom::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Natural>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_geom::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "geom" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_geom::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "geometric";
    
    return d_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the geometric distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_geom::getParameterRules(void) const 
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "The probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_geom::getTypeSpec( void ) const 
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** 
 * Set a member variable.
 * 
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_geom::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "p" ) 
    {
        p = var;
    }
    else 
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
    
}
