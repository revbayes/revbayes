#include "Transform_Vector_Invlogit.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

namespace Core = RevBayesCore;

Transform_Vector_InvLogit::Transform_Vector_InvLogit() : TypedDistribution< ModelVector<Probability> >(),
				 base_distribution( NULL )
{
    markAsTransform();
}

Transform_Vector_InvLogit::~Transform_Vector_InvLogit()
{

}

Transform_Vector_InvLogit* RevLanguage::Transform_Vector_InvLogit::clone( void ) const
{
    return new Transform_Vector_InvLogit(*this);
}

Core::TransformedVectorDistribution* RevLanguage::Transform_Vector_InvLogit::createDistribution( void ) const
{
    using namespace Transforms;

    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    std::unique_ptr<Core::TypedDistribution<Core::RbVector<double>>> vp( static_cast<Core::TypedDistribution<Core::RbVector<double>>* >( rl_vp.createDistribution() ) );

    return new Core::TransformedVectorDistribution(vp, invlogit_transform, invlogit_inverse, log_invlogit_prime);
}


/* Get Rev type of object */
const std::string& RevLanguage::Transform_Vector_InvLogit::getClassType(void)
{

    static std::string rev_type = "Transform_Vector_InvLogit";

    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Transform_Vector_InvLogit::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Probability> >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string RevLanguage::Transform_Vector_InvLogit::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "invlogit";

    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Transform_Vector_InvLogit::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	std::vector<TypeSpec> distTypes = { TypedDistribution<ModelVector<Real>>::getClassTypeSpec(), TypedDistribution<ModelVector<RealPos>>::getClassTypeSpec(), TypedDistribution<ModelVector<Probability>>::getClassTypeSpec()};

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Transform_Vector_InvLogit::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void RevLanguage::Transform_Vector_InvLogit::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< ModelVector<Probability> >::setConstParameter(name, var);
    }
}

