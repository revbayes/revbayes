/**
 * @file
 * This file contains the declaration of the uniform partitioning prior distribution, which is used create
 * random variables of UPP distributions.
 *
 * @brief Declaration of the uniform partitioning distribution
 *
 * Will Freyman 6/29/17
 *
 */

#ifndef Dist_upp_H
#define Dist_upp_H

#include "UniformPartitioningDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_upp :  public TypedDistribution< ModelVector<valType> >{
        
    public:
        Dist_upp( void );
        virtual ~Dist_upp();
        
        // Basic utility functions
        Dist_upp*                                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
//        void                                            printValue(std::ostream& o) const;                                                    //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::UniformPartitioningDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       values;
        RevPtr<const RevVariable>                       include_zero;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_upp<valType>::Dist_upp() : TypedDistribution< ModelVector<valType> >(), 
    values( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_upp<valType>::~Dist_upp()
{
    
}



template <typename valType>
RevLanguage::Dist_upp<valType>* RevLanguage::Dist_upp<valType>::clone( void ) const
{
    return new Dist_upp(*this);
}


template <typename valType>
RevBayesCore::UniformPartitioningDistribution<typename valType::valueType>* RevLanguage::Dist_upp<valType>::createDistribution( void ) const
{
	
    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* v   = static_cast<const ModelVector<valType> &>( values->getRevObject() ).getDagNode();
    bool iz = static_cast<const RlBoolean &>(include_zero->getRevObject()).getValue();
    
    RevBayesCore::UniformPartitioningDistribution<typename valType::valueType>* d		= new RevBayesCore::UniformPartitioningDistribution<typename valType::valueType>(v, iz);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_upp<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_dpp";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_upp<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_upp<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "UPP";
    
    return d_name;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_upp<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
    
        dist_member_rules.push_back( new ArgumentRule( "values", ModelVector<valType>::getClassTypeSpec(), "The possible non-zero values in the most fragmented partition of the set.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "include_zero", RlBoolean::getClassTypeSpec(), "Should the partitions also include zero values?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_upp<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_upp<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "values" )
    {
        values = var;
    }
    else if ( name == "include_zero" )
    {
        include_zero = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif
