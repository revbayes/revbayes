//
//  Dist_cppNormal.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 11/19/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_cppNormal.h"
#include "CompoundPoissonNormalDistribution.h"
#include "Natural.h"
#include "Probability.h"
#include "StochasticNode.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_cppNormal::Dist_cppNormal() : ContinuousDistribution()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_cppNormal* Dist_cppNormal::clone( void ) const
{
    
    return new Dist_cppNormal(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::CompoundPoissonNormalDistribution* Dist_cppNormal::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* rate    = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* m       = static_cast<const Real &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s       = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
    
    RevBayesCore::CompoundPoissonNormalDistribution* d        = new RevBayesCore::CompoundPoissonNormalDistribution( rate, m, s );
    
    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_cppNormal::getClassType(void)
{
    
    static std::string rev_type = "Dist_cppNormal";
    
	return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_cppNormal::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Real>::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_cppNormal::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "cppNormal";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Poisson distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_cppNormal::getParameterRules(void) const
{
    
    static MemberRules distPoisMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distPoisMemberRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The rate of the Poisson distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distPoisMemberRules.push_back( new ArgumentRule( "mu", Real::getClassTypeSpec(), "The mean of the normal distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distPoisMemberRules.push_back( new ArgumentRule( "sigma", RealPos::getClassTypeSpec(), "The standard deviation of the normal distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distPoisMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_cppNormal::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_cppNormal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else
    {
        TypedDistribution<Real>::setConstParameter(name, var);
    }
    
}
