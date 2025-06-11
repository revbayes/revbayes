#ifndef Dist_IID_H
#define Dist_IID_H

#include "IidDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_IID : public TypedDistribution< ModelVector<valType> > {
        
    public:
        Dist_IID( void );
        virtual ~Dist_IID();
        
        // Basic utility functions
        Dist_IID*                                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::IidDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       num_samples;
        RevPtr<const RevVariable>                       value_distribution;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_IID<valType>::Dist_IID() : TypedDistribution< ModelVector<valType> >(),
    num_samples( NULL ),
    value_distribution( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_IID<valType>::~Dist_IID()
{
    
}



template <typename valType>
RevLanguage::Dist_IID<valType>* RevLanguage::Dist_IID<valType>::clone( void ) const
{
    return new Dist_IID(*this);
}


template <typename valType>
RevBayesCore::IidDistribution<typename valType::valueType>* RevLanguage::Dist_IID<valType>::createDistribution( void ) const
{
    
    // get the parameters
    const Distribution& rl_vp                                           = static_cast<const Distribution &>( value_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<typename valType::valueType>* vp    = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rl_vp.createDistribution() );

    std::int64_t n = static_cast<const Natural & >( num_samples->getRevObject() ).getValue();
    
    RevBayesCore::IidDistribution<typename valType::valueType>* d = new RevBayesCore::IidDistribution<typename valType::valueType>(n, vp);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_IID<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_IID";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_IID<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector< valType > >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
template <typename valType>
std::string RevLanguage::Dist_IID<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "IID";
    
    return d_name;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_IID<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "numValues", Natural::getClassTypeSpec(), "The number of IID samples.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "valueDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The base distribution for the per category values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_IID<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_IID<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "numValues" )
    {
        num_samples = var;
    }
    else if ( name == "valueDistribution" )
    {
        value_distribution = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif

