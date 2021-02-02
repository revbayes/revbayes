//
//  Dist_outgroupBirthDeath.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 5/5/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantRateOutgroupBirthDeathProcess.h"
#include "Dist_outgroupBirthDeath.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBirthDeathProcess.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_outgroupBirthDeath::Dist_outgroupBirthDeath() : BirthDeathProcess()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_outgroupBirthDeath* Dist_outgroupBirthDeath::clone( void ) const
{
    return new Dist_outgroupBirthDeath(*this);
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
RevBayesCore::ConstantRateOutgroupBirthDeathProcess* Dist_outgroupBirthDeath::createDistribution( void ) const
{
    
    // get the parameters
    
    // the start age
    RevBayesCore::TypedDagNode<double>* ra      = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();
    // speciation rate
    RevBayesCore::TypedDagNode<double>* s       = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    // extinction rate
    RevBayesCore::TypedDagNode<double>* e       = static_cast<const RealPos &>( mu->getRevObject() ).getDagNode();
    // sampling probability (ingroup)
    RevBayesCore::TypedDagNode<double>* ri      = static_cast<const Probability &>( rhoIngroup->getRevObject() ).getDagNode();
    // sampling probability (outgroup)
    RevBayesCore::TypedDagNode<double>* ro      = static_cast<const Probability &>( rhoOutgroup->getRevObject() ).getDagNode();
    // sampling strategy (ingroup)
    const std::string& ssi                      = static_cast<const RlString &>( samplingStrategyIngroup->getRevObject() ).getValue();
    // strategy (outgroup)
    const std::string& sso                      = static_cast<const RlString &>( samplingStrategyOutgroup->getRevObject() ).getValue();
    
    // condition
    const std::string& cdt                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t         = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    std::vector<RevBayesCore::Taxon> t_out     = static_cast<const ModelVector<Taxon> &>( taxaOutgroup->getRevObject() ).getValue();
    
    // create the internal distribution object
    RevBayesCore::ConstantRateOutgroupBirthDeathProcess* d;

    d = new RevBayesCore::ConstantRateOutgroupBirthDeathProcess(ra, s, e, ri, ro, ssi, sso, cdt, t, t_out);
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_outgroupBirthDeath::getClassType( void )
{
    
    static std::string rev_type = "Dist_outgroupBirthDeath";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_outgroupBirthDeath::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( BirthDeathProcess::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_outgroupBirthDeath::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "OutgroupBDP" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_outgroupBirthDeath::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "OutgroupBirthDeath";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the constant-rate birth-death process are:
 * (1) the speciation rate lambda which must be a positive real.
 * (2) the extinction rate mu that must be a positive real.
 * (3) all member rules specified by BirthDeathProcess.
 *
 * \return The member rules.
 */
const MemberRules& Dist_outgroupBirthDeath::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "lambda",          RealPos::getClassTypeSpec(), "The constant speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu",              RealPos::getClassTypeSpec(), "The constant extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "rhoIngroup",             Probability::getClassTypeSpec(), "The ingroup sampling probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "rhoOutgroup",     Probability::getClassTypeSpec(), "The outgroup sampling probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "taxaOutgroup",    ModelVector<Taxon>::getClassTypeSpec(), "The taxa used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        std::vector<std::string> optionsStrategy;
        optionsStrategy.push_back( "uniform" );
        optionsStrategy.push_back( "diversified" );
        dist_member_rules.push_back( new OptionRule( "samplingStrategyOutgroup", new RlString("uniform"), optionsStrategy, "The sampling strategy of including outgroup taxa at the present." ) );
        dist_member_rules.push_back( new OptionRule( "samplingStrategyIngroup", new RlString("uniform"), optionsStrategy, "The sampling strategy of including ingroup taxa at the present." ) );
//        dist_member_rules.push_back( new ArgumentRule( "backbone",        TimeTree::getClassTypeSpec(), "The backbone topological constraints.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = BirthDeathProcess::getParameterRules();
        dist_member_rules.insert(dist_member_rules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_outgroupBirthDeath::getTypeSpec( void ) const
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
void Dist_outgroupBirthDeath::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "rhoIngroup" )
    {
        rhoIngroup = var;
    }
    else if ( name == "rhoOutgroup" )
    {
        rhoOutgroup = var;
    }
    else if ( name == "taxaOutgroup" )
    {
        taxaOutgroup = var;
    }
    else if ( name == "samplingStrategyOutgroup" )
    {
        samplingStrategyOutgroup = var;
    }
    else if ( name == "samplingStrategyIngroup" )
    {
        samplingStrategyIngroup = var;
    }
    else
    {
        BirthDeathProcess::setConstParameter(name, var);
    }
    
}
