#include "Transform_Logit.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

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

RevBayesCore::TransformedDistribution* Transform_Logit::createDistribution( void ) const
{
    using namespace Transforms;

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

