//
//  Func_normalize.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 9/8/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_normalize.h"
#include "ModelVector.h"
#include "NormalizeVectorFunction.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Func_normalize::Func_normalize() : TypedFunction< ModelVector<RealPos> >()
{
    
}

/* Clone object */
Func_normalize* Func_normalize::clone( void ) const {
    
    return new Func_normalize( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_normalize::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> > *params = static_cast< ModelVector<RealPos> & >( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode< double > *sum = static_cast< RealPos& >( args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::NormalizeVectorFunction *func = new RevBayesCore::NormalizeVectorFunction( params, sum );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_normalize::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<RealPos>::getClassTypeSpec(), "The vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sum", RealPos::getClassTypeSpec(), "The sum the vector will have after normalization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_normalize::getClassType(void)
{
    
    static std::string rev_type = "Func_normalize";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_normalize::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_normalize::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "normalize";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_normalize::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
