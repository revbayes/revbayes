//
//  Func__or.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/16/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func__or.h"

#include "ArgumentRule.h"
#include "LogicalOrFunction.h"
#include "RbUtil.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Func__or::Func__or() : TypedFunction<RlBoolean>() {
    
}

/* Clone object */
Func__or* Func__or::clone( void ) const
{
    
    return new Func__or( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Boolean>* Func__or::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode<RevBayesCore::Boolean>* leftVal = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<RevBayesCore::Boolean>* rightVal = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::LogicalOrFunction *func = new RevBayesCore::LogicalOrFunction( leftVal, rightVal );

    return func;
}


/** Get argument rules */
const ArgumentRules& Func__or::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "", RlBoolean::getClassTypeSpec(), "first value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "", RlBoolean::getClassTypeSpec(), "second value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func__or::getClassType(void)
{
    
    static std::string rev_type = "Func__or";
    
	return rev_type; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func__or::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "or";
    
    return f_name;
}


/** Get class type spec describing type of object */
const TypeSpec& Func__or::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& Func__or::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

