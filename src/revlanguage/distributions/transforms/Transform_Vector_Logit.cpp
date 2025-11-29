#include "Transform_Vector_Logit.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace Core = RevBayesCore;

Transform_Vector_Logit::Transform_Vector_Logit() : TypedDistribution< ModelVector<Real> >(),
				 base_distribution( NULL )
{
    markAsTransform();
}

Transform_Vector_Logit::~Transform_Vector_Logit()
{

}

Transform_Vector_Logit* RevLanguage::Transform_Vector_Logit::clone( void ) const
{
    return new Transform_Vector_Logit(*this);
}

std::optional<double> logit_transform(double x);
std::optional<double> logit_inverse(double x);
std::optional<double> log_logit_prime(double x);

Core::TransformedVectorDistribution* RevLanguage::Transform_Vector_Logit::createDistribution( void ) const
{
    // get the parameters
    const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
    std::unique_ptr<Core::TypedDistribution<Core::RbVector<double>>> vp( static_cast<Core::TypedDistribution<Core::RbVector<double>>* >( rl_vp.createDistribution() ) );

    return new Core::TransformedVectorDistribution(vp, logit_transform, logit_inverse, log_logit_prime);
}


/* Get Rev type of object */
const std::string& RevLanguage::Transform_Vector_Logit::getClassType(void)
{

    static std::string rev_type = "Transform_Vector_Logit";

    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Transform_Vector_Logit::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Real> >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string RevLanguage::Transform_Vector_Logit::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "logit";

    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Transform_Vector_Logit::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
	// Should we add 
	std::vector<TypeSpec> distTypes = { TypedDistribution<ModelVector<Probability>>::getClassTypeSpec(), TypedDistribution<Simplex>::getClassTypeSpec() };

        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Transform_Vector_Logit::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}



/** Set a member variable */
void RevLanguage::Transform_Vector_Logit::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< ModelVector<Real> >::setConstParameter(name, var);
    }
}

