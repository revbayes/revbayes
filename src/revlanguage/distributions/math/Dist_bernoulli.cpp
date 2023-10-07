#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BernoulliDistribution.h"
#include "Dist_bernoulli.h"
#include "Probability.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "Natural.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;

Dist_bernoulli::Dist_bernoulli() : TypedDistribution<Natural>()
{
    
    setGuiDistributionName("Bernoulli");
    setGuiDistributionToolTip("Bernoulli distribution for 0 or 1 random variables");
}


Dist_bernoulli::~Dist_bernoulli()
{
    
}



Dist_bernoulli* Dist_bernoulli::clone( void ) const
{
    
    return new Dist_bernoulli(*this);
}


RevBayesCore::BernoulliDistribution* Dist_bernoulli::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* q   = static_cast<const Probability &>( p->getRevObject() ).getDagNode();
    RevBayesCore::BernoulliDistribution* d  = new RevBayesCore::BernoulliDistribution( q );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_bernoulli::getClassType(void)
{
    
    static std::string rev_type = "Dist_bernoulli";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_bernoulli::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_bernoulli::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "bernoulli";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_bernoulli::getParameterRules(void) const
{
    
    static MemberRules distBernMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        distBernMemberRules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "The probability of success.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distBernMemberRules;
}


const TypeSpec& Dist_bernoulli::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_bernoulli::printValue(std::ostream& o) const
{
    
    o << "Bernoulli(p=";
    if ( p != NULL )
    {
        o << p->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_bernoulli::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if ( name == "p" ) 
    {
        p = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
}
