#ifndef Dist_WeightedSample_H
#define Dist_WeightedSample_H

#include "WeightedSampleDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_WeightedSample :  public TypedDistribution< ModelVector<valType> >{
        
    public:
        Dist_WeightedSample( void );
        virtual ~Dist_WeightedSample();
        
        // Basic utility functions
        Dist_WeightedSample*                            clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::WeightedSampleDistribution<typename valType::valueType>*  createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       baseDistribution;
        RevPtr<const RevVariable>                       numSamples;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RealPos.h"
#include "RbVector.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_WeightedSample<valType>::Dist_WeightedSample() : TypedDistribution< ModelVector<valType> >(),
baseDistribution( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_WeightedSample<valType>::~Dist_WeightedSample()
{
    
}



template <typename valType>
RevLanguage::Dist_WeightedSample<valType>* RevLanguage::Dist_WeightedSample<valType>::clone( void ) const
{
    return new Dist_WeightedSample(*this);
}


template <typename valType>
RevBayesCore::WeightedSampleDistribution<typename valType::valueType>* RevLanguage::Dist_WeightedSample<valType>::createDistribution( void ) const
{
    
    RevBayesCore::WeightedSampleDistribution<typename valType::valueType>* d = NULL;
    RevBayesCore::TypedDistribution<typename valType::valueType>* g0  = NULL;
    size_t ns = 0;
    std::vector< RevBayesCore::TypedDistribution<typename valType::valueType>* > g0_vector;
    
    // get the parameters
    if (baseDistribution->getRevObject().isType( WorkspaceVector<TypedDistribution<valType> >::getClassTypeSpec() ) ) {

        // get the workspace vector
        const RevObject& obj = baseDistribution->getRevObject();
        const WorkspaceVector<TypedDistribution<valType> >& obj_vector = static_cast<const WorkspaceVector<TypedDistribution<valType> >& >( obj );
        ns = obj_vector.size();
        
        // build our vector of typed distributions
        g0_vector = std::vector< RevBayesCore::TypedDistribution<typename valType::valueType>* >( ns, NULL );
        for (size_t i = 0; i < ns; i++)
        {
            const Distribution& rlDistribution = static_cast<const Distribution &>( obj_vector[i] );
            g0_vector[i] = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );
        }
    }
    else {
        ns                                 = size_t( static_cast<const Natural &>( numSamples->getRevObject() ).getValue() );
        const Distribution& rlDistribution = static_cast<const Distribution &>( baseDistribution->getRevObject() );
        g0                                 = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );
//        d = new RevBayesCore::WeightedSampleDistribution<typename valType::valueType>(g0, ns);
    }
    
    // create the weighted same distribution
    d = new RevBayesCore::WeightedSampleDistribution<typename valType::valueType>(g0, g0_vector, ns);
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_WeightedSample<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_WeightedSample";
    
    return rev_type;
}


/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_WeightedSample<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_WeightedSample<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "WeightedSample";
    
    return d_name;
}



template <typename valType>
RevLanguage::MethodTable RevLanguage::Dist_WeightedSample<valType>::getDistributionMethods( void ) const
{
    // MJL: Note, getDistributionMethods() is currently disabled! Need to fix this.
    MethodTable methods;
    return methods;
    
//    const Distribution& rlDistribution = static_cast<const Distribution &>( baseDistribution->getRevObject() );
//    return rlDistribution.getDistributionMethods();
}

/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_WeightedSample<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_empirical_sample_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        std::vector<TypeSpec> baseDistributionTypes;
        baseDistributionTypes.push_back( TypedDistribution<valType>::getClassTypeSpec() );
        baseDistributionTypes.push_back( WorkspaceVector<TypedDistribution<valType> >::getClassTypeSpec() );
//        dist_member_rules.push_back( new ArgumentRule( "Q", rateMatrixTypes, "The global, branch-specific or site-mixture rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        dist_empirical_sample_member_rules.push_back( new ArgumentRule( "baseDistribution", baseDistributionTypes, "The base distribution(s) for the sampled values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_empirical_sample_member_rules.push_back( new ArgumentRule( "numSamples", Natural::getClassTypeSpec(), "The number of samples in the weighted average.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural() ) );
        
        rules_set = true;
    }
    
    return dist_empirical_sample_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_WeightedSample<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_WeightedSample<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "baseDistribution" )
    {
        baseDistribution = var;
    }
    else if ( name == "numSamples" )
    {
        numSamples = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif
