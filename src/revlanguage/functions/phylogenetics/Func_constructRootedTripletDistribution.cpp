//
//  Func_constructRootedTripletDistribution.cpp
//  RevBayesCore
//
//  Created by Bastien Boussau on 8/7/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func_constructRootedTripletDistribution.h"

#include <cstddef>

#include "ModelVector.h"
#include "RlBranchLengthTree.h"
#include "RlDeterministicNode.h"
#include "RlRootedTripletDistribution.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTaxon.h"
#include "RootedTripletDistributionFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "RlTree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** default constructor */
Func_constructRootedTripletDistribution::Func_constructRootedTripletDistribution( void ) : TypedFunction<RootedTripletDistribution>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_constructRootedTripletDistribution* Func_constructRootedTripletDistribution::clone( void ) const {
    
    return new Func_constructRootedTripletDistribution( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RootedTripletDistribution >* Func_constructRootedTripletDistribution::createFunction( void ) const
{
    
    RevBayesCore::RootedTripletDistributionFunction* f = NULL;
    
    
    if ( this->args[1].getVariable()->getRequiredTypeSpec().isDerivedOf( ModelVector< RlString >::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector< std::string > >* sn = static_cast<const ModelVector< RlString > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        f = new RevBayesCore::RootedTripletDistributionFunction( sn );
    }
    else if ( this->args[1].getVariable()->getRequiredTypeSpec().isDerivedOf( ModelVector< Taxon >::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector< RevBayesCore::Taxon > >* t = static_cast<const ModelVector< Taxon > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        f = new RevBayesCore::RootedTripletDistributionFunction( t );
    }

    if ( this->args[0].getVariable()->getRequiredTypeSpec().isDerivedOf( ModelVector< TimeTree >::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector< RevBayesCore::Tree > >* gTrees = static_cast<const ModelVector< TimeTree > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        f->setTrees(gTrees);
    }
    else if ( this->args[0].getVariable()->getRequiredTypeSpec().isDerivedOf( ModelVector< BranchLengthTree >::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector< RevBayesCore::Tree > >* gTrees = static_cast<const ModelVector< BranchLengthTree > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        f->setTrees(gTrees);
    }

    if ( this->args[2].getVariable()->getRequiredTypeSpec().isDerivedOf( RlBoolean::getClassTypeSpec() ) )
    {
        
        // Sebastian: Currently unused
//        RevBayesCore::TypedDagNode< RevBayesCore::Boolean >* t = static_cast<const RlBoolean &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
        // TODO: Bastien, this isn't working (Sebastian)
//        f->setRecordBranchLengths( t );
    }

    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_constructRootedTripletDistribution::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        std::vector<TypeSpec> treeTypes;
        treeTypes.push_back( ModelVector< TimeTree >::getClassTypeSpec() );
        treeTypes.push_back( ModelVector< BranchLengthTree >::getClassTypeSpec() );


        argumentRules.push_back( new ArgumentRule( "geneTrees",         Tree::getClassTypeSpec(),                       "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "speciesNames",      ModelVector< RlString >::getClassTypeSpec(),    "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "keepBranchLengths", RlBoolean::getClassTypeSpec(),                  "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_constructRootedTripletDistribution::getClassType(void)
{
    
    static std::string rev_type = "Func_constructRootedTripletDistribution";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_constructRootedTripletDistribution::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_constructRootedTripletDistribution::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "rootedTripletDist";
    
    return f_name;
}


const TypeSpec& Func_constructRootedTripletDistribution::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
