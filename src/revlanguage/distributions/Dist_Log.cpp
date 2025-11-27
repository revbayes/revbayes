#include "Dist_Log.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

Dist_Log::Dist_Log() : TypedDistribution< RealPos >(),
    log_distribution( NULL )
{
    
}

Dist_Log::~Dist_Log()
{
    
}

Dist_Log* RevLanguage::Dist_Log::clone( void ) const
{
    return new Dist_Log(*this);
}

RevBayesCore::TransformedDistribution* RevLanguage::Dist_Log::createDistribution( void ) const
{
    using namespace Transforms;

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( log_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TransformedDistribution* d = new RevBayesCore::TransformedDistribution(*vp, exp_transform, exp_inverse, log_exp_prime);

    delete vp;
    
    return d;
}


/* Get Rev type of object */
const std::string& RevLanguage::Dist_Log::getClassType(void)
{
    
    static std::string rev_type = "Dist_Log";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Dist_Log::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< RealPos >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string RevLanguage::Dist_Log::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Log";
    
    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Dist_Log::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "logDistribution", distTypes, "The distribution in log-space.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Dist_Log::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
void RevLanguage::Dist_Log::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "logDistribution" )
    {
        log_distribution = var;
    }
    else
    {
        TypedDistribution< RealPos >::setConstParameter(name, var);
    }
}

