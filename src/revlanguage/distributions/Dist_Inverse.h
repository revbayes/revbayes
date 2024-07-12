#ifndef Dist_Inverse_H
#define Dist_Inverse_H

// Need to include headers for all possible values of valType
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "Simplex.h"
#include "IidDistribution.h"
#include "ModelVector.h"
#include "TypedDistribution.h"
#include "InverseDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {

    template<typename valType>
    class Dist_Inverse : public TypedDistribution< valType > {
        
    public:
        Dist_Inverse( void );
        virtual ~Dist_Inverse();
        
        // Basic utility functions
        Dist_Inverse*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::InverseDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       dist;
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
using namespace RevLanguage;

template<typename valType>
Dist_Inverse<valType>::Dist_Inverse<valType>() : TypedDistribution< valType >(),
    dist( NULL )
{
    
}

Dist_Inverse::~Dist_Inverse()
{
    
}

Dist_Inverse* RevLanguage::Dist_Inverse::clone( void ) const
{
    return new Dist_Inverse(*this);
}


template<typename valType>
RevBayesCore::InverseDistribution* RevLanguage::Dist_Inverse::createDistribution( void ) const
{
    
    // get the parameters
    const Distribution& rl_vp = static_cast<const Distribution &>( dist->getRevObject() );
    
    // Cast the single distribution to the specific type expected by InverseDistribution
    RevBayesCore::TypedDistribution<typename valType::valueType>* vp  = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rl_vp.createDistribution() );

    // Create an instance of InverseDistribution using the single distribution
    RevBayesCore::InverseDistribution<typename valType::valueType>* d = new RevBayesCore::InverseDistribution<typename valType::valueType>(vp);

    delete vp;
    
    return d;

}


/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_Inverse<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_Inverse";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_Inverse<valType>::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
template <typename valType>
std::vector<std::string> Dist_Inverse<valType>::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "inv" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
template <typename valType>
std::string RevLanguage::Dist_Inverse<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "inverse";
    
    return d_name;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_Inverse<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "distribution",
         TypedDistribution<valType>::getClassTypeSpec(), "The distribution to invert.",
          ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const RevLanguage::TypeSpec& RevLanguage::Dist_Inverse::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_Inverse<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "distribution" )
    {
        inverse_distribution = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif // Dist_Inverse_H
