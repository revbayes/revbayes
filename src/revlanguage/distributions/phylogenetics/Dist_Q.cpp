#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BetaDistribution.h"
#include "Dist_Q.h"
#include "RlSimplex.h"
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

Dist_Q::Dist_Q() : TypedDistribution<RateGenerator>()
{
    
    setGuiDistributionName("Q");
    setGuiDistributionToolTip("Q distribution for random variables on rate matrices");
}


Dist_Q::~Dist_Q()
{
    
}



Dist_Q* Dist_Q::clone( void ) const
{
    
    return new Dist_Q(*this);
}


RevBayesCore::QDistribution* Dist_Q::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* a   = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* b   = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    RevBayesCore::QDistribution* d          = new RevBayesCore::QDistribution(a, b);
    
    return d;
}



RateMatrix* Dist_Q::createRandomVariable(void) const
{
    
    RevBayesCore::TypedDistribution<RevBayesCore::RateGenerator>* d = createDistribution();
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rv  = new StochasticNode("", d, this->clone() );
    
    return new RateMatrix(rv);
}



/* Get Rev type of object */
const std::string& Dist_Q::getClassType(void)
{
    
    static std::string rev_type = "Dist_Q";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_Q::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Probability>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_Q::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "beta";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_Q::getParameterRules(void) const
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


const TypeSpec& Dist_Q::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_Q::printValue(std::ostream& o) const
{
    
    o << "Q(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ", prob=";
    if ( rho != NULL )
    {
        o << rho->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_Q::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else
    {
        TypedDistribution<RateGenerator>::setConstParameter(name, var);
    }
}
