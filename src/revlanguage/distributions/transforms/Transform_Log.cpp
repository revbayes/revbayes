#include "Transform_Log.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Log::Transform_Log() : TypedDistribution< Real >(),
                                 base_distribution( NULL )
{
    markAsTransform();
}

Transform_Log::~Transform_Log()
{
}

Transform_Log* Transform_Log::clone( void ) const
{
    return new Transform_Log(*this);
}

std::optional<double> log_transform(double x)
{
    if (x > 0)
        return log(x);
    else
        return {}; // out of range
}

std::optional<double> log_inverse(double x)
{
    return exp(x);
}

std::optional<double> log_log_prime(double x)
{
    // y = log(x)
    // dy/dx = 1/x
    // log(dy/dx) = -log(x)

    if (x > 0)
        return -log(x);
    else
        return {}; // out of range
}


RevBayesCore::TransformedDistribution* Transform_Log::createDistribution( void ) const
{

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TransformedDistribution* d = new RevBayesCore::TransformedDistribution(*vp, log_transform, log_inverse, log_log_prime);

    delete vp;

    return d;
}


/* Get Rev type of object */
const std::string& Transform_Log::getClassType(void)
{

    static std::string rev_type = "Transform_Log";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Log::getClassTypeSpec(void)
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
std::string Transform_Log::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "log";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Log::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        std::vector<TypeSpec> distTypes = { TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec() };

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Transform_Log::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Log::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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

