//
//  Func__unot.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/16/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func__unot.h"

#include "ArgumentRule.h"
#include "LogicalUnotFunction.h"
#include "RbUtil.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Func__unot::Func__unot() : TypedFunction<RlBoolean>() {
    
}

/* Clone object */
Func__unot* Func__unot::clone( void ) const {
    
    return new Func__unot( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Boolean>* Func__unot::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode<RevBayesCore::Boolean>* val = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::LogicalUnotFunction *func = new RevBayesCore::LogicalUnotFunction( val );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func__unot::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) {
        
        argumentRules.push_back( new ArgumentRule( "", RlBoolean::getClassTypeSpec(), "The expression.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func__unot::getClassType(void)
{
    
    static std::string rev_type = "Func__unot";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func__unot::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func__unot::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "unot";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func__unot::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

