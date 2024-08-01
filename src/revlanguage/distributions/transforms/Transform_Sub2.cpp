#include "Transform_Sub2.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Sub2::Transform_Sub2() : TypedDistribution< Real >()
{
    markAsTransform();
}

Transform_Sub2::~Transform_Sub2()
{
}

Transform_Sub2* Transform_Sub2::clone( void ) const
{
    return new Transform_Sub2(*this);
}

RevBayesCore::TransformedDistribution* Transform_Sub2::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* fir        = static_cast<const Real &>( first->getRevObject() ).getDagNode();

    const Distribution& rl_vp                      = static_cast<const Distribution &>( second_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    // compute the transforms
    auto sub2_transform = [=](double x) -> optional<double> { return fir->getValue() - x; };
    auto sub2_inverse = sub2_transform;
    auto log_sub2_prime = [=](double x) -> optional<double> { return 0; };

    RevBayesCore::TransformedDistribution* dist = new RevBayesCore::TransformedDistribution(*vp, sub2_transform, sub2_inverse, log_sub2_prime, {fir});

    delete vp;

    return dist;
}


/* Get Rev type of object */
const std::string& Transform_Sub2::getClassType(void)
{

    static std::string rev_type = "Transform_Sub2";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Sub2::getClassTypeSpec(void)
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
std::string Transform_Sub2::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "_sub";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Sub2::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "first", Real::getClassTypeSpec()   , "The amount sub2ed to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "secondDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Transform_Sub2::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Sub2::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "secondDistribution" )
    {
        second_distribution = var;
    }
    else if ( name == "first" )
    {
        first = var;
    }
    else
    {
        TypedDistribution< Real >::setConstParameter(name, var);
    }
}

