#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_exponential.h"
#include "ExponentialDistribution.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "RlPositiveContinuousDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_exponential::Dist_exponential() : PositiveContinuousDistribution() {
    
}


Dist_exponential::~Dist_exponential() {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_exponential* Dist_exponential::clone( void ) const {
    return new Dist_exponential(*this);
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
RevBayesCore::ExponentialDistribution* Dist_exponential::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<double>* l     = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* o     = static_cast<const Real    &>( offset->getRevObject() ).getDagNode();
    RevBayesCore::ExponentialDistribution* d  = new RevBayesCore::ExponentialDistribution( l, o );
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_exponential::getClassType(void)
{
    
    static std::string rev_type = "Dist_exponential";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_exponential::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_exponential::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "exp" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_exponential::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "exponential";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_exponential::getParameterRules(void) const
{
    
    static MemberRules distExpMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        distExpMemberRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate parameter ( rate = 1/mean ).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        distExpMemberRules.push_back( new ArgumentRule( "offset", RealPos::getClassTypeSpec(), "The offset of the distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        
        rules_set = true;
    }
    
    return distExpMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_exponential::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_exponential::printValue(std::ostream& o) const {
    
    o << " exponential(lambda=";
    if ( lambda != NULL ) {
        o << lambda->getName();
    } else {
        o << "?";
    }
    o << ", offset=";
    if ( offset != NULL ) {
        o << offset->getName();
    } else {
        o << "?";
    }
    o << ")";
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
void Dist_exponential::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "offset" )
    {
        offset = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
