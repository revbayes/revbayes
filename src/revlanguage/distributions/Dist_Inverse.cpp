#include "Dist_Inverse.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
using namespace RevLanguage;

template<typename valType>
Dist_Inverse<valType>::Dist_Inverse<valType>() : TypedDistribution< valType >(),
    inverse_distribution( NULL )
{
    
}

Dist_Inverse::~Dist_Inverse()
{
    
}

Dist_Inverse* RevLanguage::Dist_Inverse::clone( void ) const
{
    return new Dist_Inverse(*this);
}


RevBayesCore::InverseDistribution* RevLanguage::Dist_Inverse::createDistribution( void ) const
{
    
    // get the parameters
    const Distribution& rl_vp                          = static_cast<const Distribution &>( inverse_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<valType>* vp  = static_cast<RevBayesCore::TypedDistribution<valType>* >( rl_vp.createDistribution() );

    RevBayesCore::InverseDistribution* d = new RevBayesCore::InverseDistribution(*vp);

    delete vp;
    
    return d;
}


/* Get Rev type of object */
const std::string& RevLanguage::Dist_Inverse::getClassType(void)
{
    
    static std::string rev_type = "Dist_Inverse";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::Dist_Inverse::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_Inverse::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Inverse";
    
    return d_name;
}


/** Return member rules (no members) */
const RevLanguage::MemberRules& RevLanguage::Dist_Inverse::getParameterRules(void) const
{
    static MemberRules dist_member_rules;    
    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Dist_Inverse::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
void RevLanguage::Dist_Inverse::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "inverseDistribution" )
    {
        inverse_distribution = var;
    }
    else
    {
        TypedDistribution< valType >::setConstParameter(name, var);
    }
}
