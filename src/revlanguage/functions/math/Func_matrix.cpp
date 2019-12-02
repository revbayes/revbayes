//
//  Func_matrix.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 7/5/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_matrix.h"

#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "ModelVector.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "RlSimplex.h"
#include "MatrixFunction.h"


using namespace RevLanguage;


/** Default constructor */
Func_matrix::Func_matrix( void ) : TypedFunction< MatrixReal >()
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_matrix* Func_matrix::clone( void ) const
{
    return new Func_matrix( *this );
}


/** Execute function: Construct simplex from RealPos values. */
RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_matrix::createFunction( void ) const
{
//    std::vector<std::vector<const RevBayesCore::TypedDagNode<double>* > > params;
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* mtx;
    mtx = static_cast< const ModelVector<ModelVector<Real> >& >( args[0].getVariable()->getRevObject() ).getDagNode();
//
//    for ( size_t i = 0; i < args[.size(); i++ )
//    {
//        params.push_back( std::vector<const RevBayesCore::TypedDagNode<double>* >() );
//        for ( size_t j = 0; j < args[i].size(); j++ )
//        {
//            const Real &val = static_cast<const Real &>( args[i][j].getVariable()->getRevObject() );
//        }
//        params.push_back( val.getDagNode() );
//    }

    RevBayesCore::MatrixFunction* func = new RevBayesCore::MatrixFunction( mtx );
    
//    RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* func = NULL;
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_matrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<ModelVector<Real> >::getClassTypeSpec(), "x", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static) */
const std::string& Func_matrix::getClassType( void )
{
    static std::string rev_type = "Func_matrix";
    
    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_matrix::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_matrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "matrix";
    
    return f_name;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_matrix::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

