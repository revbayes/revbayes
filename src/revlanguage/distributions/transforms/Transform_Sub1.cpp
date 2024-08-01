#include "Transform_Sub1.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Sub1::Transform_Sub1() : TypedDistribution< Real >()
{
    markAsTransform();
}

Transform_Sub1::~Transform_Sub1()
{
}

Transform_Sub1* Transform_Sub1::clone( void ) const
{
    return new Transform_Sub1(*this);
}

RevBayesCore::TransformedDistribution* Transform_Sub1::createDistribution( void ) const
{
    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( first_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TypedDagNode<double>* sec        = static_cast<const Real &>( second->getRevObject() ).getDagNode();

    auto sub1_transform = [=](double x) -> optional<double> { return x - sec->getValue(); };
    auto sub1_inverse = [=](double x) -> optional<double> { return x + sec->getValue(); };
    auto log_sub1_prime = [=](double x) -> optional<double> { return 0; };

    RevBayesCore::TransformedDistribution* dist = new RevBayesCore::TransformedDistribution(*vp, sub1_transform, sub1_inverse, log_sub1_prime, {sec});

    delete vp;

    return dist;
}


/* Get Rev type of object */
const std::string& Transform_Sub1::getClassType(void)
{

    static std::string rev_type = "Transform_Sub1";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Sub1::getClassTypeSpec(void)
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
std::string Transform_Sub1::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "_sub";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Sub1::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "firstDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "second", Real::getClassTypeSpec()   , "The amount sub1ed to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Transform_Sub1::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Sub1::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "firstDistribution" )
    {
        first_distribution = var;
    }
    else if ( name == "second" )
    {
        second = var;
    }
    else
    {
        TypedDistribution< Real >::setConstParameter(name, var);
    }
}

