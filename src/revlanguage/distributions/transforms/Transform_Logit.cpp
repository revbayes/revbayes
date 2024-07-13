#include "Transform_Logit.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Logit::Transform_Logit() : TypedDistribution< Real >(),
                                 base_distribution( NULL )
{
    markAsTransform();
}

Transform_Logit::~Transform_Logit()
{
}

Transform_Logit* Transform_Logit::clone( void ) const
{
    return new Transform_Logit(*this);
}

// This can handle x=exp(-1000), but not x=1 - exp(-1000)

std::optional<double> logit_transform(double x)
{
    if (x > 0 and x < 1)
    {
	// log( x / (1-x) )
	return log(x) - log1p(-x);
    }
    else
        return {}; // out of range
}

std::optional<double> logit_inverse(double x)
{
    // Two forms of the same expression to avoid overflow with exp(larg number)
    if (x < 0)
	return exp(x)/(1+exp(x));
    else
	return 1/(exp(-x)+1);
}

std::optional<double> log_logit_prime(double x)
{
    // y = log(x/(1-x))
    // dy/dx = 1/(x/(1-x)) * ((1-x)(1) - x(-1))/((1-x)**2)
    //       = (1-x)/x     * 1/((1-x)**2)
    //       = 1/(x * (1-x) )
    // logit(dy/dx) = -log(x) - log(1-x)

    if (x > 0 and x < 1)
        return -log(x) - log1p(-x);
    else
        return {}; // out of range
}


RevBayesCore::TransformedDistribution* Transform_Logit::createDistribution( void ) const
{

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TransformedDistribution* d = new RevBayesCore::TransformedDistribution(*vp, logit_transform, logit_inverse, log_logit_prime);

    delete vp;

    return d;
}


/* Get Rev type of object */
const std::string& Transform_Logit::getClassType(void)
{

    static std::string rev_type = "Transform_Logit";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Logit::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< Real >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Transform_Logit::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "logit";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Logit::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        std::vector<TypeSpec> distTypes = { TypedDistribution<Probability>::getClassTypeSpec() };

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Transform_Logit::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Logit::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< Real >::setConstParameter(name, var);
    }
}

