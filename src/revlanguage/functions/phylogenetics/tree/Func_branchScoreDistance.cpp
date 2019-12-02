//
//  Func_branchScoreDistance.cpp
//  RevBayesCore
//
//  Created by Bastien Boussau on 8/7/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func_branchScoreDistance.h"
#include "ModelVector.h"
#include "RlTimeTree.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "BranchScoreDistanceStatistic.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_branchScoreDistance::Func_branchScoreDistance( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_branchScoreDistance* Func_branchScoreDistance::clone( void ) const {
    
    return new Func_branchScoreDistance( *this );
}


RevBayesCore::TypedFunction< double >* Func_branchScoreDistance::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau2 = static_cast<const Tree&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::BranchScoreDistanceStatistic* f = new RevBayesCore::BranchScoreDistanceStatistic( tau, tau2 );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_branchScoreDistance::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree1", Tree::getClassTypeSpec(), "The first tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "tree2", Tree::getClassTypeSpec(), "The second tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_branchScoreDistance::getClassType(void)
{
    
    static std::string rev_type = "Func_branchScoreDistance";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_branchScoreDistance::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_branchScoreDistance::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "branchScoreDistance";
    
    return f_name;
}


const TypeSpec& Func_branchScoreDistance::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
