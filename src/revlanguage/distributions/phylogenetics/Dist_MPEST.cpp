//
//  Dist_MPEST.cpp
//  RevBayes
//
//  Created by Bastien on 30/10/14.
//  Copyright (c) 2014 Bastien Boussau. All rights reserved.
//


#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_MPEST.h"
#include "MPEST.h"
#include "RlTimeTree.h"
#include "RlRootedTripletDistribution.h"
#include "RootedTripletDistribution.h"
#include "StochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_MPEST::Dist_MPEST() : TypedDistribution<RootedTripletDistribution>()
{
    
}

/**
 * Clone the object
 *
 * \return a clone of the object.
 */

Dist_MPEST* Dist_MPEST::clone(void) const {
    
    return new Dist_MPEST(*this);
    
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
RevBayesCore::MPEST* Dist_MPEST::createDistribution( void ) const
{
    
    // Get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::RootedTripletDistribution>* st = static_cast<const RootedTripletDistribution &>( speciesTree->getRevObject() ).getDagNode();
    //RevBayesCore::TypedDagNode<RevBayesCore::RootedTripletDistribution>* gt = static_cast<const RootedTripletDistribution &>( geneTrees->getRevObject() ).getDagNode();

    const bool useSp = static_cast<const RlBoolean& >( useSpecies->getRevObject() ).getValue();
    
    RevBayesCore::MPEST*   d = new RevBayesCore::MPEST( st, useSp );
    
    //We can't do that because we can't simulate using MPEST    d->redrawValue();
    
    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_MPEST::getClassType(void)
{
    
    static std::string rev_type = "Dist_MPEST";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_MPEST::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_MPEST::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MPEST";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Multispecies Coalescent process are:
 * (1) Species tree rooted triplet distribution.
 * (2) Whether we use species names or taxons.
 *
 * \return The member rules.
 */
const MemberRules& Dist_MPEST::getMemberRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree", RootedTripletDistribution::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        //MPESTMemberRules.push_back( new ArgumentRule( "geneTrees", RootedTripletDistribution::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        memberRules.push_back( new ArgumentRule( "useSpecies", RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
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
void Dist_MPEST::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "speciesTree" )
    {
        speciesTree = var;
    }
  /*  else if ( name == "geneTrees" )
    {
        geneTrees = var;
    }*/
    else if ( name == "useSpecies" )
    {
        useSpecies = var;
    }
    else
    {
        TypedDistribution<RootedTripletDistribution>::setConstParameter(name, var);
    }
    
}


/** Get type spec */
const TypeSpec& Dist_MPEST::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



