//
//  Dist_varianceGamma.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/8/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_varianceGamma.h"
#include "VarianceGammaDistribution.h"
#include "Real.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_varianceGamma::Dist_varianceGamma() : ContinuousDistribution()
{
    
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
RevBayesCore::VarianceGammaDistribution* Dist_varianceGamma::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* m        = static_cast<const Real &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* k        = static_cast<const RealPos &>( kappa->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* t        = static_cast<const RealPos &>( tau->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ti       = static_cast<const RealPos &>( tau->getRevObject() ).getDagNode();
    RevBayesCore::VarianceGammaDistribution*   d = new RevBayesCore::VarianceGammaDistribution(m, k, t, ti);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_varianceGamma* Dist_varianceGamma::clone( void ) const
{
    
    return new Dist_varianceGamma(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_varianceGamma::getClassType(void)
{
    
    static std::string rev_type = "Dist_varianceGamma";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_varianceGamma::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_varianceGamma::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "VG" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_varianceGamma::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "VarianceGamma";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the mean of the distribution.
 * (2) the standard deviation.
 *
 * \return The member rules.
 */
const MemberRules& Dist_varianceGamma::getParameterRules(void) const
{
    
    static MemberRules distVarianceGammaMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distVarianceGammaMemberRules.push_back( new ArgumentRule( "mu", Real::getClassTypeSpec()   , "The mean parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        distVarianceGammaMemberRules.push_back( new ArgumentRule( "kappa"  , RealPos::getClassTypeSpec(), "The standard deviation parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        distVarianceGammaMemberRules.push_back( new ArgumentRule( "tau"  , RealPos::getClassTypeSpec(), "The standard deviation parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        distVarianceGammaMemberRules.push_back( new ArgumentRule( "time"  , RealPos::getClassTypeSpec(), "The duration of time for the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );


        rules_set = true;
    }
    
    return distVarianceGammaMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_varianceGamma::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_varianceGamma::printValue(std::ostream& o) const
{
    
    o << " varianceGamma(mu=";
    if ( mu != NULL )
    {
        o << mu->getName();
    }
    else
    {
        o << "?";
    }
    o << ", kappa=";
    if ( kappa != NULL )
    {
        o << kappa->getName();
    }
    else
    {
        o << "?";
    }
    o << ", tau=";
    if ( tau != NULL )
    {
        o << tau->getName();
    }
    else
    {
        o << "?";
    }
    o << ", time=";
    if ( time != NULL )
    {
        o << time->getName();
    }
    else
    {
        o << "?";
    }

    o << ")";
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
void Dist_varianceGamma::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "kappa" )
    {
        kappa = var;
    }
    else if ( name == "tau" )
    {
        tau = var;
    }
    else if ( name == "time" )
    {
        time = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
