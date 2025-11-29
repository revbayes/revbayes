#include "Transform_Exp.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

Transform_Exp::Transform_Exp() : TypedDistribution< RealPos >(),
				 base_distribution( NULL )
{
    markAsTransform();
}

Transform_Exp::~Transform_Exp()
{
    
}

Transform_Exp* RevLanguage::Transform_Exp::clone( void ) const
{
    return new Transform_Exp(*this);
}

std::optional<double> exp_transform(double x)
{
    return exp(x);
}

std::optional<double> exp_inverse(double x)
{
    if (x > 0)
	return log(x);
    else
	return {}; // out of range
}

std::optional<double> log_exp_prime(double x)
{
    // y = exp(x)
    // dy/dx = exp(x)
    // log(dy/dx) = x;

    return x;
}


RevBayesCore::TransformedDistribution* RevLanguage::Transform_Exp::createDistribution( void ) const
{
    
    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TransformedDistribution* d = new RevBayesCore::TransformedDistribution(*vp, exp_transform, exp_inverse, log_exp_prime);

    delete vp;
    
    return d;
}


/* Get Rev type of object */
const std::string& RevLanguage::Transform_Exp::getClassType(void)
{
    
    static std::string rev_type = "Transform_Exp";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Transform_Exp::getClassTypeSpec(void)
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
std::string RevLanguage::Transform_Exp::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "exp";
    
    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Transform_Exp::getParameterRules(void) const
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


const RevLanguage::TypeSpec& RevLanguage::Transform_Exp::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
void RevLanguage::Transform_Exp::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< RealPos >::setConstParameter(name, var);
    }
}

