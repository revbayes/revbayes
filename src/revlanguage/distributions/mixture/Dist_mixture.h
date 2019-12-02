#ifndef Dist_mixture_H
#define Dist_mixture_H

#include "MixtureDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
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
     * @since 2014-08-14, version 1.0
     *
     */
    template <typename valType>
    class Dist_mixture : public TypedDistribution< valType >{
        
    public:
                                                        Dist_mixture( void );
        virtual                                        ~Dist_mixture();
        
        // Basic utility functions
        Dist_mixture*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)

        
        // Distribution functions you have to override
        RevBayesCore::MixtureDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       values;
        RevPtr<const RevVariable>                       probabilities;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_mixture<valType>::Dist_mixture() : TypedDistribution< valType >(),
    values( NULL ),
    probabilities( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_mixture<valType>::~Dist_mixture()
{
    
}



template <typename valType>
RevLanguage::Dist_mixture<valType>* RevLanguage::Dist_mixture<valType>::clone( void ) const
{
    
    return new Dist_mixture(*this);
}


template <typename valType>
RevBayesCore::MixtureDistribution<typename valType::valueType>* RevLanguage::Dist_mixture<valType>::createDistribution( void ) const
{
	
    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* v   = static_cast<const ModelVector<valType> &>( values->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*                      p   = static_cast<const Simplex &>( probabilities->getRevObject() ).getDagNode();
    
    RevBayesCore::MixtureDistribution<typename valType::valueType>* d		= new RevBayesCore::MixtureDistribution<typename valType::valueType>(v,p);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_mixture<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_mixture";
    
	return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixture<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_mixture<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Mixture";
    
    return d_name;
}

template <typename valType>
RevLanguage::MethodTable RevLanguage::Dist_mixture<valType>::getDistributionMethods( void ) const
{

    MethodTable methods = TypedDistribution<valType>::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_allocation_index_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_mixture, Natural >( "getAllocationIndex", this->variable, get_allocation_index_arg_rules, true ) );

    return methods;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_mixture<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "values"       , ModelVector<valType>::getClassTypeSpec(), "The potential values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "probabilities", Simplex::getClassTypeSpec()             , "The probabilitoes for each value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_mixture<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_mixture<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "values" )
    {
        values = var;
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
