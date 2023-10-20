//
//  Func_heterotachyIndex.cpp
//  RevBayesCore
//
//  Created by April Wright, Laura Mulvey and Basanta Khakurel on 10/20/23.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func_HeterotachyIndex.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "HeterotachyIndex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** default constructor */
Func_heterotachyIndex::Func_heterotachyIndex( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_heterotachyIndex* Func_heterotachyIndex::clone( void ) const {
    
    return new Func_heterotachyIndex( *this );
}


RevBayesCore::TypedFunction< double >* Func_heterotachyIndex::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau1 = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau2 = static_cast<const Tree&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::heterotachyIndex* f = new RevBayesCore::heterotachyIndex( tau1, tau2 );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_heterotachyIndex::getArgumentRules( void ) const
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


const std::string& Func_heterotachyIndex::getClassType(void)
{
    
    static std::string rev_type = "Func_heterotachyIndex";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_heterotachyIndex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_heterotachyIndex::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "heterotachyIndex";
    
    return f_name;
}


const TypeSpec& Func_heterotachyIndex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
