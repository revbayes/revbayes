#ifndef Dist_EmpiricalSample_H
#define Dist_EmpiricalSample_H

#include "EmpiricalSampleDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_EmpiricalSample :  public TypedDistribution< ModelVector<valType> >{
        
    public:
        Dist_EmpiricalSample( void );
        virtual ~Dist_EmpiricalSample();
        
        // Basic utility functions
        Dist_EmpiricalSample*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::EmpiricalSampleDistribution<typename valType::valueType>*  createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:

        RevPtr<const RevVariable>                       baseDistribution;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlDistributionMemberFunction.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_EmpiricalSample<valType>::Dist_EmpiricalSample() : TypedDistribution< ModelVector<valType> >(),
baseDistribution( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_EmpiricalSample<valType>::~Dist_EmpiricalSample()
{
    
}



template <typename valType>
RevLanguage::Dist_EmpiricalSample<valType>* RevLanguage::Dist_EmpiricalSample<valType>::clone( void ) const
{
    return new Dist_EmpiricalSample(*this);
}


template <typename valType>
RevBayesCore::EmpiricalSampleDistribution<typename valType::valueType>* RevLanguage::Dist_EmpiricalSample<valType>::createDistribution( void ) const
{
    
    // get the parameters
    const Distribution& rlDistribution                                          = static_cast<const Distribution &>( baseDistribution->getRevObject() );
    RevBayesCore::TypedDistribution<typename valType::valueType>* g0            = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );

    RevBayesCore::EmpiricalSampleDistribution<typename valType::valueType>* d   = new RevBayesCore::EmpiricalSampleDistribution<typename valType::valueType>(g0);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_EmpiricalSample<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_EmpiricalSample";
    
    return rev_type;
}


/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_EmpiricalSample<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_EmpiricalSample<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "EmpiricalSample";
    
    return d_name;
}



template <typename valType>
RevLanguage::MethodTable RevLanguage::Dist_EmpiricalSample<valType>::getDistributionMethods( void ) const
{
    
    const Distribution& rlDistribution = static_cast<const Distribution &>( baseDistribution->getRevObject() );
    
    
    MethodTable methods = rlDistribution.getDistributionMethods();
    
    // member functions
    ArgumentRules* sample_prob_arg_rules = new ArgumentRules();
    sample_prob_arg_rules->push_back( new ArgumentRule( "log",       RlBoolean::getClassTypeSpec(), "If we should return the log-transformed probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );
    sample_prob_arg_rules->push_back( new ArgumentRule( "normalize", RlBoolean::getClassTypeSpec(), "If we should normalize the probabilities.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );
    methods.addFunction( new DistributionMemberFunction<Dist_EmpiricalSample<valType>, ModelVector<Real> >( "getSampleProbabilities", this->variable, sample_prob_arg_rules, true   ) );
    
    return methods;
}

/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_EmpiricalSample<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_empirical_sample_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_empirical_sample_member_rules.push_back( new ArgumentRule( "baseDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The base distribution for the sampled values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_empirical_sample_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_EmpiricalSample<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_EmpiricalSample<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "baseDistribution" )
    {
        baseDistribution = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif
