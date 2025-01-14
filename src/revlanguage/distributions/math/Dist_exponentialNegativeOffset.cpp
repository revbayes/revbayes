#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_exponentialNegativeOffset.h"
#include "ExponentialDistribution.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Dist_exponentialNegativeOffset::Dist_exponentialNegativeOffset() : ContinuousDistribution() {
    
}


Dist_exponentialNegativeOffset::~Dist_exponentialNegativeOffset() {
    
}



Dist_exponentialNegativeOffset* Dist_exponentialNegativeOffset::clone( void ) const {
    return new Dist_exponentialNegativeOffset(*this);
}


RevBayesCore::ExponentialDistribution* Dist_exponentialNegativeOffset::createDistribution( void ) const {
    // get the parameters
    RevBayesCore::TypedDagNode<double>* l     = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* o     = static_cast<const Real    &>( offset->getRevObject() ).getDagNode();
    RevBayesCore::ExponentialDistribution* d  = new RevBayesCore::ExponentialDistribution( l, o );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_exponentialNegativeOffset::getClassType(void)
{
    
    static std::string rev_type = "Dist_exponentialNegativeOffset";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_exponentialNegativeOffset::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_exponentialNegativeOffset::getDistributionFunctionAliases( void ) const
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
std::string Dist_exponentialNegativeOffset::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "exponential";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_exponentialNegativeOffset::getParameterRules(void) const
{
    
    static MemberRules distExpMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distExpMemberRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate parameter ( rate = 1/mean ) of the distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distExpMemberRules.push_back( new ArgumentRule( "offset", Real::getClassTypeSpec()   , "The offset parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distExpMemberRules;
}


const TypeSpec& Dist_exponentialNegativeOffset::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_exponentialNegativeOffset::printValue(std::ostream& o) const
{
    
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


/** Set a member variable */
void Dist_exponentialNegativeOffset::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
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
