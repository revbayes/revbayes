#include "Dist_MultivariateLog.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

namespace Core = RevBayesCore;

Dist_MultivariateLog::Dist_MultivariateLog() : TypedDistribution< ModelVector<RealPos> >(),
    log_distribution( NULL )
{
    
}

Dist_MultivariateLog::~Dist_MultivariateLog()
{
    
}

Dist_MultivariateLog* RevLanguage::Dist_MultivariateLog::clone( void ) const
{
    return new Dist_MultivariateLog(*this);
}

Core::TransformedVectorDistribution* RevLanguage::Dist_MultivariateLog::createDistribution( void ) const
{
    using namespace Transforms;
    using Core::RbVector;

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( log_distribution->getRevObject() );
    std::unique_ptr<Core::TypedDistribution<Core::RbVector<double>>> vp( static_cast<Core::TypedDistribution<Core::RbVector<double>>* >( rl_vp.createDistribution() ) );;

    return new Core::TransformedVectorDistribution( vp,
						    exp_transform,
						    exp_inverse,
						    log_exp_prime);
}


/* Get Rev type of object */
const std::string& RevLanguage::Dist_MultivariateLog::getClassType(void)
{
    
    static std::string rev_type = "Dist_MultivariateLog";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Dist_MultivariateLog::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<RealPos> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string RevLanguage::Dist_MultivariateLog::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Log";
    
    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Dist_MultivariateLog::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<ModelVector<Real>>::getClassTypeSpec(), TypedDistribution<ModelVector<RealPos>>::getClassTypeSpec(), TypedDistribution<ModelVector<Probability>>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "logDistribution", distTypes, "The distribution in log-space.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Dist_MultivariateLog::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
void RevLanguage::Dist_MultivariateLog::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "logDistribution" )
    {
        log_distribution = var;
    }
    else
    {
        TypedDistribution< ModelVector<RealPos> >::setConstParameter(name, var);
    }
}

