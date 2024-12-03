#include "Transform_Shift.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevLanguage;

Transform_Shift::Transform_Shift() : TypedDistribution< Real >(),
                                     base_distribution( NULL )
{
    markAsTransform();
}

Transform_Shift::~Transform_Shift()
{
}

Transform_Shift* Transform_Shift::clone( void ) const
{
    return new Transform_Shift(*this);
}

RevBayesCore::TransformedDistribution* Transform_Shift::createDistribution( void ) const
{
    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

    RevBayesCore::TypedDagNode<double>* d           = static_cast<const Real &>( delta->getRevObject() ).getDagNode();

    auto shift_transform = [=](double x) -> optional<double> { return x + d->getValue(); };
    auto shift_inverse = [=](double x) -> optional<double> { return x - d->getValue(); };
    auto log_shift_prime = [=](double x) -> optional<double> { return 0; };

    RevBayesCore::TransformedDistribution* dist = new RevBayesCore::TransformedDistribution(*vp, shift_transform, shift_inverse, log_shift_prime, {d});

    delete vp;

    return dist;
}


/* Get Rev type of object */
const std::string& Transform_Shift::getClassType(void)
{

    static std::string rev_type = "Transform_Shift";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Transform_Shift::getClassTypeSpec(void)
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
std::string Transform_Shift::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "shift";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Transform_Shift::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<Real>::getClassTypeSpec(), TypedDistribution<RealPos>::getClassTypeSpec(), TypedDistribution<Probability>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "delta", Real::getClassTypeSpec()   , "The amount added to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Transform_Shift::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void Transform_Shift::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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

