//
//  Func_vectorFlatten.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//


#include "ArgumentRule.h"
#include "Func_vectorFlatten.h"
#include "Real.h"
#include "VectorFlattenFunction.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTypedFunction.h"


using namespace RevLanguage;


/** Default constructor */
Func_vectorFlatten::Func_vectorFlatten( void ) :
TypedFunction<ModelVector<Real> >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_vectorFlatten* Func_vectorFlatten::clone( void ) const
{
    return new Func_vectorFlatten( *this );
}


/** Execute function: Compute simplex from vector of RealPos values */
RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_vectorFlatten::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* vec;
    vec = static_cast< const ModelVector<ModelVector<Real> >& >( args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::VectorFlattenFunction* func = new RevBayesCore::VectorFlattenFunction( vec );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_vectorFlatten::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<ModelVector<Real> >::getClassTypeSpec(), "A vector of a vector.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static version) */
const std::string& Func_vectorFlatten::getClassType( void )
{
    static std::string rev_type = "Func_vectorFlatten";
    
	return rev_type;
}


/** Get Rev type spec of object (static version) */
const TypeSpec& Func_vectorFlatten::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_vectorFlatten::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "vectorFlatten";
    
    return f_name;
}


/** Get Rev type spec from an instance of the object */
const TypeSpec& Func_vectorFlatten::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

