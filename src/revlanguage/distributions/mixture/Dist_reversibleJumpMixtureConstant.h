#ifndef Dist_reversibleJumpMixtureConstant_H
#define Dist_reversibleJumpMixtureConstant_H

#include "ReversibleJumpMixtureConstantDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the ReversibleJumpMixtureConstantDistribution
     *
     * The RevLanguage wrapper of the mixture distribution connects
     * the variables/parameters of the distribution and creates the internal ReversibleJumpMixtureConstantDistribution object.
     * Please read the ReversibleJumpMixtureConstantDistribution.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-14, version 1.0
     *
     */
    template <typename valType>
    class Dist_reversibleJumpMixtureConstant :  public TypedDistribution< valType >{
        
    public:
        Dist_reversibleJumpMixtureConstant( void );
        virtual                                        ~Dist_reversibleJumpMixtureConstant();
        
        // Basic utility functions
        Dist_reversibleJumpMixtureConstant*             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::ReversibleJumpMixtureConstantDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       constValue;
        RevPtr<const RevVariable>                       baseDistribution;
        RevPtr<const RevVariable>                       probability;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Probability.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::Dist_reversibleJumpMixtureConstant() : TypedDistribution< valType >(),
    constValue( NULL ),
    baseDistribution( NULL ),
    probability( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::~Dist_reversibleJumpMixtureConstant()
{
    
}



template <typename valType>
RevLanguage::Dist_reversibleJumpMixtureConstant<valType>* RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::clone( void ) const
{
    
    return new Dist_reversibleJumpMixtureConstant(*this);
}


template <typename valType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<typename valType::valueType>* RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode< typename valType::valueType >*      cv              = static_cast<const valType &>( constValue->getRevObject() ).getDagNode();
    const Distribution&                                             rlDistribution  = static_cast<const Distribution &>( baseDistribution->getRevObject() );
    RevBayesCore::TypedDistribution<typename valType::valueType>*   g0              = dynamic_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );

    // error checking
    if ( g0 == NULL )
    {
        throw RbException("Could not cast distribution in reversible jump mixture.");
    }
    RevBayesCore::TypedDagNode< double >*                           p               = static_cast<const Probability &>( probability->getRevObject() ).getDagNode();
    
    RevBayesCore::ReversibleJumpMixtureConstantDistribution<typename valType::valueType>* d		= new RevBayesCore::ReversibleJumpMixtureConstantDistribution<typename valType::valueType>(cv,g0,p);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_reversibleJumpMixtureConstant";
    
    return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
template <typename valType>
std::vector<std::string> RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "RJMixture" );
    
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
std::string RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "ReversibleJumpMixture";
    
    return d_name;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        //dist_member_rules.push_back( new ArgumentRule( "constantValue", ModelVector<valType>::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        dist_member_rules.push_back( new ArgumentRule( "constantValue", valType::getClassTypeSpec(), "The fixed value this distribution can take.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "baseDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The distribution from which the value is alternatively drawn.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "p", Probability::getClassTypeSpec(), "The probability of being the fixed value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_reversibleJumpMixtureConstant<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "constantValue" )
    {
        constValue = var;
    }
    else if ( name == "baseDistribution" )
    {
        baseDistribution = var;
    }
    else if ( name == "p" )
    {
        probability = var;
    }
    else
    {
        TypedDistribution< valType >::setConstParameter(name, var);
    }
}


#endif
