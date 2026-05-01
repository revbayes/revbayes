#ifndef Dist_mixtureVector_H
#define Dist_mixtureVector_H

#include "VectorMixtureDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "ModelVector.h"
#include "TypeSpec.h"

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
    class Dist_mixtureVector : public TypedDistribution< ModelVector<valType> >{
        
    public:
        Dist_mixtureVector( void );
        virtual                                        ~Dist_mixtureVector();
        
        // Basic utility functions
        Dist_mixtureVector*                             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::VectorMixtureDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       base_distributions;
        RevPtr<const RevVariable>                       probabilities;
        RevPtr<const RevVariable>                       num_values;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_mixtureVector<valType>::Dist_mixtureVector() : TypedDistribution< ModelVector<valType> >(),
    base_distributions( NULL ),
    probabilities( NULL ),
    num_values( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_mixtureVector<valType>::~Dist_mixtureVector()
{
    
}



template <typename valType>
RevLanguage::Dist_mixtureVector<valType>* RevLanguage::Dist_mixtureVector<valType>::clone( void ) const
{
    
    return new Dist_mixtureVector(*this);
}


template <typename valType>
RevBayesCore::VectorMixtureDistribution< typename valType::valueType >* RevLanguage::Dist_mixtureVector<valType>::createDistribution( void ) const
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
    long n = static_cast<const Natural &>( num_values->getRevObject() ).getValue();

    RevBayesCore::VectorMixtureDistribution<typename valType::valueType>* d = new RevBayesCore::VectorMixtureDistribution<typename valType::valueType>(bd,p,n);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_mixtureVector<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_mixtureVector";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixtureVector<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_mixtureVector<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MixtureVector";
    
    return d_name;
}

template <typename valType>
RevLanguage::MethodTable RevLanguage::Dist_mixtureVector<valType>::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution< ModelVector<valType> >::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_mixture_prob_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_mixtureVector, Simplex >( "getMixtureProbabilities", this->variable, get_mixture_prob_arg_rules, true ) );
    
    return methods;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_mixtureVector<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "baseDistributions", WorkspaceVector< TypedDistribution<valType> >::getClassTypeSpec(), "The base distribution for the per category values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "probabilities"    , Simplex::getClassTypeSpec(), "The probabilities for each value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "numValues",        Natural::getClassTypeSpec(), "The number of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixtureVector<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_mixtureVector<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "baseDistributions" )
    {
        base_distributions = var;
    }
    else if ( name == "probabilities" )
    {
        probabilities = var;
    }
    else if ( name == "numValues" )
    {
        num_values = var;
    }
    else
    {
        TypedDistribution< ModelVector<valType> >::setConstParameter(name, var);
    }
}


#endif
