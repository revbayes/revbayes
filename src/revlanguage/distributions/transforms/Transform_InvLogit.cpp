#include "Transform_InvLogit.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_InvLogit::Transform_InvLogit() : TypedDistribution< Probability >(),
                                 base_distribution( NULL )
{
    markAsTransform();
}

Transform_InvLogit::~Transform_InvLogit()
{
}

Transform_InvLogit* Transform_InvLogit::clone( void ) const
{
    return new Transform_InvLogit(*this);
}

std::optional<double> invlogit_transform(double x)
{
    // Two forms of the same expression to avoid overflow with exp(larg number)
    if (x < 0)
	return exp(x)/(1+exp(x));
    else
	return 1/(exp(-x)+1);
}

std::optional<double> invlogit_inverse(double x)
{
    if (x > 0 and x < 1)
    {
	// log( x / (1-x) )
	return log(x) - log1p(-x);
    }
    else
        return {}; // out of range
}

std::optional<double> log_invlogit_prime(double x)
{
    if (x > 0)
    {
	// y = 1/(exp(-x) + 1)
	// dy/dx = -1/((exp(-x)+1)^2) * exp(-x) * -1
	//       = exp(-x) / (1 + exp(-x))^2
	// log(dy/dx) = -x - 2*log(1+exp(-x))

        return - x - 2*log1p(exp(-x));
    }
    else
    {
	// y = exp(x)/(exp(x) + 1)
	// dy/dx = exp(x) / (1+exp(x))^2
	// log(dy/dx) = x  - 2log1p(exp(x))
    
	return   x - 2*log1p(exp(x));
    }
}


RevBayesCore::TransformedDistribution* Transform_InvLogit::createDistribution( void ) const
{

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TransformedDistribution* d = new RevBayesCore::TransformedDistribution(*vp, invlogit_transform, invlogit_inverse, log_invlogit_prime);

    delete vp;

    return d;
}


/* Get Rev type of object */
const std::string& Transform_InvLogit::getClassType(void)
{

    static std::string rev_type = "Transform_InvLogit";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_InvLogit::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< Probability >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Transform_InvLogit::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "invlogit";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_InvLogit::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Transform_InvLogit::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_InvLogit::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< Probability >::setConstParameter(name, var);
    }
}

