#include "Transform_Mul.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

using namespace RevLanguage;
namespace Core = RevBayesCore;

Transform_Mul::Transform_Mul() : TypedDistribution< Real >(),
                                     base_distribution( NULL )
{
    markAsTransform();
}

Transform_Mul::~Transform_Mul()
{
}

Transform_Mul* Transform_Mul::clone( void ) const
{
    return new Transform_Mul(*this);
}

Core::TransformedDistribution* Transform_Mul::createDistribution( void ) const
{
    using namespace Transforms;

    // get the parameters
    const Distribution& rl_vp               = static_cast<const Distribution &>( base_distribution->getRevObject() );
    Core::TypedDistribution<double>* vp     = static_cast<Core::TypedDistribution<double>* >( rl_vp.createDistribution() );

    Core::TypedDagNode<double>* l           = static_cast<const Real &>( delta->getRevObject() ).getDagNode();

    Core::TransformedDistribution* dist = new Core::TransformedDistribution(*vp, mul_transform, mul_inverse, log_mul_prime, {l});

    delete vp;

    return dist;
}


/* Get Rev type of object */
const std::string& Transform_Mul::getClassType(void)
{

    static std::string rev_type = "Transform_Mul";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Mul::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec ( TypedDistribution< Real >::getClassTypeSpec() );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Transform_Mul::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "_mul";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Mul::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "delta", Real::getClassTypeSpec()   , "The amount muled to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Transform_Mul::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Mul::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
    }
    else
    {
        TypedDistribution< Real >::setConstParameter(name, var);
    }
}

