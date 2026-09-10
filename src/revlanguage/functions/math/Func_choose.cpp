//
//  Func_choose.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/21/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "Func_choose.h"

#include "ChooseFunction.h"
#include "Natural.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_choose::Func_choose( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_choose* Func_choose::clone( void ) const
{
    
    return new Func_choose( *this );
}


RevBayesCore::TypedFunction<std::int64_t>* Func_choose::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<std::int64_t>* n = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* k = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::ChooseFunction* f = new RevBayesCore::ChooseFunction( n, k );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_choose::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "n"    , Natural::getClassTypeSpec(), "Number of total elements.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "k"    , Natural::getClassTypeSpec(), "Number of chosen elements.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_choose::getClassType(void)
{
    
    static std::string rev_type = "Func_choose";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_choose::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_choose::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "choose";
    
    return f_name;
}


const TypeSpec& Func_choose::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
