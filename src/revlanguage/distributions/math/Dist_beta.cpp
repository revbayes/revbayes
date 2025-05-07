#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BetaDistribution.h"
#include "Dist_beta.h"
#include "RealPos.h"
#include "Probability.h"
#include "RlContinuousStochasticNode.h"
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
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;

Dist_beta::Dist_beta() : ProbabilityContinuousDistribution()
{
    
    setGuiDistributionName("Beta");
    setGuiDistributionToolTip("Beta distribution for random variables on the interval [0,1]");
}


Dist_beta::~Dist_beta()
{
    
}



Dist_beta* Dist_beta::clone( void ) const
{
    
    return new Dist_beta(*this);
}


RevBayesCore::BetaDistribution* Dist_beta::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* a   = static_cast<const RealPos &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* b   = static_cast<const RealPos &>( beta->getRevObject() ).getDagNode();
    RevBayesCore::BetaDistribution* d       = new RevBayesCore::BetaDistribution(a, b);
    
    return d;
}



Probability* Dist_beta::createRandomVariable(void) const
{
    
    RevBayesCore::ContinuousDistribution* d = createDistribution();
    RevBayesCore::TypedDagNode<double>* rv  = new ContinuousStochasticNode("", d, this->clone() );
    
    return new Probability(rv);
}



/* Get Rev type of object */
const std::string& Dist_beta::getClassType(void)
{
    
    static std::string rev_type = "Dist_beta";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_beta::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ProbabilityContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_beta::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "beta";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_beta::getParameterRules(void) const
{
    
    static MemberRules distUnifMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distUnifMemberRules.push_back( new ArgumentRule( "alpha", RealPos::getClassTypeSpec(), "The alpha shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distUnifMemberRules.push_back( new ArgumentRule( "beta" , RealPos::getClassTypeSpec(), "The beta shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distUnifMemberRules;
}


const TypeSpec& Dist_beta::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_beta::printValue(std::ostream& o) const
{
    
    o << "beta(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ", beta=";
    if ( beta != NULL )
    {
        o << beta->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_beta::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "beta" )
    {
        beta = var;
    }
    else
    {
        ProbabilityContinuousDistribution::setConstParameter(name, var);
    }
}
