#include "Transform_Scale_Pos.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Scale_Pos::Transform_Scale_Pos() : TypedDistribution< RealPos >(),
                                     base_distribution( NULL )
{
    markAsTransform();
}

Transform_Scale_Pos::~Transform_Scale_Pos()
{
}

Transform_Scale_Pos* Transform_Scale_Pos::clone( void ) const
{
    return new Transform_Scale_Pos(*this);
}

RevBayesCore::TransformedDistribution* Transform_Scale_Pos::createDistribution( void ) const
{
    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TypedDagNode<double>* l           = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();

    auto scale_transform = [=](double x) -> optional<double> { return x * l->getValue(); };
    auto scale_inverse = [=](double x) -> optional<double> { return x / l->getValue(); };
    auto log_scale_prime = [=](double x) -> optional<double> { return l->getValue(); };

    RevBayesCore::TransformedDistribution* dist = new RevBayesCore::TransformedDistribution(*vp, scale_transform, scale_inverse, log_scale_prime, {l});

    delete vp;

    return dist;
}


/* Get Rev type of object */
const std::string& Transform_Scale_Pos::getClassType(void)
{

    static std::string rev_type = "Transform_Scale_Pos";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Scale_Pos::getClassTypeSpec(void)
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
std::string Transform_Scale_Pos::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "scale";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Scale_Pos::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()   , "The factor multiplied by the base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Transform_Scale_Pos::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Scale_Pos::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else
    {
        TypedDistribution< RealPos >::setConstParameter(name, var);
    }
}

