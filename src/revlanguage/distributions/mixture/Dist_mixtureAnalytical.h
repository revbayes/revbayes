#ifndef Dist_mixtureAnalytical_H
#define Dist_mixtureAnalytical_H

#include "AnalyticalMixtureDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "ModelVector.h"
#include "TypeSpec.h"

#include <vector>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the MixtureDistribution
     *
     * The RevLanguage wrapper of the mixture distribution connects
     * the variables/parameters of the distribution and creates the internal MixtureDistribution object.
     * Please read the MixtureDistribution.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2019-08-28, version 1.0.12
     *
     */
    template <typename valType>
    class Dist_mixtureAnalytical : public TypedDistribution< valType >{
        
    public:
        Dist_mixtureAnalytical( void );
        virtual                                        ~Dist_mixtureAnalytical();
        
        // Basic utility functions
        Dist_mixtureAnalytical*                         clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::AnalyticalMixtureDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       base_distributions;
        RevPtr<const RevVariable>                       probabilities;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_mixtureAnalytical<valType>::Dist_mixtureAnalytical() : TypedDistribution< valType >(),
    base_distributions( NULL ),
    probabilities( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_mixtureAnalytical<valType>::~Dist_mixtureAnalytical()
{
    
}



template <typename valType>
RevLanguage::Dist_mixtureAnalytical<valType>* RevLanguage::Dist_mixtureAnalytical<valType>::clone( void ) const
{
    
    return new Dist_mixtureAnalytical(*this);
}


template <typename valType>
RevBayesCore::AnalyticalMixtureDistribution< typename valType::valueType >* RevLanguage::Dist_mixtureAnalytical<valType>::createDistribution( void ) const
{
    
    // get the parameters
    const WorkspaceVector<Distribution>& rl_bd = static_cast<const WorkspaceVector<Distribution > &>( base_distributions->getRevObject() );
    std::vector<RevBayesCore::TypedDistribution<typename valType::valueType>* > bd;
    for (size_t i=0; i<rl_bd.size();++i)
    {
        RevBayesCore::TypedDistribution<typename valType::valueType>* tmp = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rl_bd[i].createDistribution() );
        bd.push_back( tmp );
    }
    
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* p = static_cast<const Simplex &>( probabilities->getRevObject() ).getDagNode();
    
    RevBayesCore::AnalyticalMixtureDistribution<typename valType::valueType>* d = new RevBayesCore::AnalyticalMixtureDistribution<typename valType::valueType>(bd,p);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_mixtureAnalytical<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_mixtureAnalytical";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixtureAnalytical<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
    
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
std::string RevLanguage::Dist_mixtureAnalytical<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MixtureAnalytical";
    
    return d_name;
}

template <typename valType>
RevLanguage::MethodTable RevLanguage::Dist_mixtureAnalytical<valType>::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution< valType >::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_mixture_prob_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_mixtureAnalytical, Simplex >( "getMixtureProbabilities", this->variable, get_mixture_prob_arg_rules, true ) );
    
    return methods;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_mixtureAnalytical<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "baseDistributions", WorkspaceVector< TypedDistribution<valType> >::getClassTypeSpec(), "The base distribution for the per category values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "probabilities"    , Simplex::getClassTypeSpec(), "The probabilities for each value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixtureAnalytical<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_mixtureAnalytical<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "baseDistributions" )
    {
        base_distributions = var;
    }
    else if ( name == "probabilities" )
    {
        probabilities = var;
    }
    else
    {
        TypedDistribution< valType >::setConstParameter(name, var);
    }
}


#endif
