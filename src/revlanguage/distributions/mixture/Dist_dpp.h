/**
 * @file
 * This file contains the declaration of the Dirichlet process prior distribution, which is used create
 * random variables of DPP distributions.
 *
 * @brief Declaration of the Dirichlet distribution
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-04-20 04:06:14 +0200 (Fri, 20 Apr 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since Version 1.0, 2012-08-06
 *
 * $Id: Func__add.h 1406 2012-04-20 02:06:14Z hoehna $
 */

#ifndef Dist_dpp_H
#define Dist_dpp_H

#include "DirichletProcessPriorDistribution.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_dpp :  public TypedDistribution< ModelVector<valType> >{
        
    public:
        Dist_dpp( void );
        virtual ~Dist_dpp();
        
        // Basic utility functions
        Dist_dpp*                                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
//        void                                            printValue(std::ostream& o) const;                                                    //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::DirichletProcessPriorDistribution<typename valType::valueType>*            createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       concentration;
        RevPtr<const RevVariable>                       baseDistribution;
        RevPtr<const RevVariable>                       numElements;
        
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
RevLanguage::Dist_dpp<valType>::Dist_dpp() : TypedDistribution< ModelVector<valType> >(), 
    concentration( NULL ),
    baseDistribution( NULL ),
    numElements( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_dpp<valType>::~Dist_dpp()
{
    
}



template <typename valType>
RevLanguage::Dist_dpp<valType>* RevLanguage::Dist_dpp<valType>::clone( void ) const
{
    return new Dist_dpp(*this);
}


template <typename valType>
RevBayesCore::DirichletProcessPriorDistribution<typename valType::valueType>* RevLanguage::Dist_dpp<valType>::createDistribution( void ) const
{
	
    // get the parameters
    const Distribution& rlDistribution									= static_cast<const Distribution &>( baseDistribution->getRevObject() );
    RevBayesCore::TypedDistribution<typename valType::valueType>* g0    = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );    
    RevBayesCore::TypedDagNode<double>* cp								= static_cast<const RealPos &>( concentration->getRevObject() ).getDagNode();

    int ne = (int)static_cast<const Natural &>( numElements->getRevObject() ).getValue();
    RevBayesCore::DirichletProcessPriorDistribution<typename valType::valueType>* d		= new RevBayesCore::DirichletProcessPriorDistribution<typename valType::valueType>(g0, cp, ne);
    
    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_dpp<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_dpp";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_dpp<valType>::getClassTypeSpec(void)
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
std::string RevLanguage::Dist_dpp<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "DPP";
    
    return d_name;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_dpp<valType>::getParameterRules(void) const
{
    
    static MemberRules distDPPMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
    
        distDPPMemberRules.push_back( new ArgumentRule( "concentration"   , RealPos::getClassTypeSpec()                   , "The concentration parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distDPPMemberRules.push_back( new ArgumentRule( "baseDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The base distribution for the per category values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distDPPMemberRules.push_back( new ArgumentRule( "numElements"     , Natural::getClassTypeSpec()                   , "The number of elements drawn from this distribution.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distDPPMemberRules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_dpp<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_dpp<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "concentration" )
    {
        concentration = var;
    }
    else if ( name == "baseDistribution" )
    {
        baseDistribution = var;
    }
    else if ( name == "numElements" )
    {
        numElements = var;
    }
    else
    {
        TypedDistribution< ModelVector< valType > >::setConstParameter(name, var);
    }
}


#endif
