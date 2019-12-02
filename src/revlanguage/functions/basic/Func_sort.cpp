//
//  Func_sort.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 9/8/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Func_sort.h"
#include "ArgumentRule.h"
#include "RbUtil.h"
#include "RevNullObject.h"
#include "Real.h"
#include "RlBoolean.h"
#include "SortVectorFunction.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Func_sort::Func_sort() : TypedFunction< ModelVector< Real > >()
{
    
}

/* Clone object */
Func_sort* Func_sort::clone( void ) const
{
    
    return new Func_sort( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_sort::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* vec = static_cast<const ModelVector<Real> &>( args[0].getVariable()->getRevObject() ).getDagNode();
    bool ascending = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

    RevBayesCore::SortVectorFunction *func = new RevBayesCore::SortVectorFunction( vec, ascending );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_sort::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "The vector to be sorted.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "ascending", RlBoolean::getClassTypeSpec(), "Sort in ascending order?.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_sort::getClassType(void)
{
    
    static std::string rev_type = "Func_sort";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_sort::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_sort::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "sort";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_sort::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
