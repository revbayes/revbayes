//
//  Dist_sampledSpeciationBirthDeathProcess.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/10/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_sampledSpeciationBirthDeathProcess.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlDistributionMemberFunction.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "SampledSpeciationBirthDeathProcess.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;


Dist_sampledSpeciationBirthDeathProcess::Dist_sampledSpeciationBirthDeathProcess() : TypedDistribution<TimeTree>()
{
    
}


Dist_sampledSpeciationBirthDeathProcess::~Dist_sampledSpeciationBirthDeathProcess()
{
    
}



Dist_sampledSpeciationBirthDeathProcess* Dist_sampledSpeciationBirthDeathProcess::clone( void ) const
{
    
    return new Dist_sampledSpeciationBirthDeathProcess( *this );
}


RevBayesCore::SampledSpeciationBirthDeathProcess* Dist_sampledSpeciationBirthDeathProcess::createDistribution( void ) const
{
    
    // Get the parameters
    RevBayesCore::TypedDagNode<double>* ra   = static_cast<const RealPos &>( root_age->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* sp   = static_cast<const RealPos &>( speciation->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ex   = static_cast<const RealPos &>( extinction->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rh   = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    
    RevBayesCore::SampledSpeciationBirthDeathProcess*   d = new RevBayesCore::SampledSpeciationBirthDeathProcess( ra, sp, ex, rh, t );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_sampledSpeciationBirthDeathProcess::getClassType(void)
{
    
    static std::string rev_type = "Dist_sampledSpeciationBirthDeathProcess";
    
    return rev_type;
}


/* Get class type spec describing type of object. TODO: Check if the correct parent is TypedDistribution or Distribution */
const TypeSpec& Dist_sampledSpeciationBirthDeathProcess::getClassTypeSpec(void)
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
std::string Dist_sampledSpeciationBirthDeathProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "SSBDP";
    
    return d_name;
}


MethodTable Dist_sampledSpeciationBirthDeathProcess::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();
    
    // member functions
    ArgumentRules* numEventsArgRules = new ArgumentRules();
    //    parentArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_sampledSpeciationBirthDeathProcess, ModelVector<Natural> >( "numberEvents", variable, numEventsArgRules   ) );
    methods.addFunction( new DistributionMemberFunction<Dist_sampledSpeciationBirthDeathProcess, ModelVector<ModelVector<RealPos> > >( "eventTimes", variable, numEventsArgRules   ) );

    
    return methods;
}



/* Return member rules */
const MemberRules& Dist_sampledSpeciationBirthDeathProcess::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "rootAge"   , RealPos::getClassTypeSpec()              , "The age of the root."                        , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "lambda"    , RealPos::getClassTypeSpec() , "The vector of speciation rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "mu"        , RealPos::getClassTypeSpec() , "The vector of extinction rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rho"       , Probability::getClassTypeSpec()          , "The taxon sampling probability."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "taxa"      , ModelVector<Taxon>::getClassTypeSpec()   , "The taxa used for simulation."               , ArgumentRule::BY_VALUE                , ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_sampledSpeciationBirthDeathProcess::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_sampledSpeciationBirthDeathProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "rootAge" )
    {
        root_age = var;
    }
    else if ( name == "lambda" )
    {
        speciation = var;
    }
    else if ( name == "mu" )
    {
        extinction = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}


