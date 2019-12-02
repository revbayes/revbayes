//
//  Func_empiricalQuantile.cpp
//  revbayes
//
//  Created by Lyndon Coghill on 9/22/16.
//  Copyright © 2016 Lyndon Coghill. All rights reserved.
//

#include "EmpiricalQuantileFunction.h"
#include "Func_empiricalQuantile.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Probability.h"

using namespace RevLanguage;

/** default constructor */
Func_empiricalQuantile::Func_empiricalQuantile( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_empiricalQuantile* Func_empiricalQuantile::clone( void ) const
{
    
    return new Func_empiricalQuantile( *this );
}


RevBayesCore::TypedFunction<double>* Func_empiricalQuantile::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* k = static_cast<const Probability &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::EmpiricalQuantileFunction* f = new RevBayesCore::EmpiricalQuantileFunction( arg, k );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_empiricalQuantile::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<Real>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "k", Probability::getClassTypeSpec(), "The kth quantile of vector x.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_empiricalQuantile::getClassType(void)
{
    
    static std::string rev_type = "Func_empiricalQuantile";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_empiricalQuantile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_empiricalQuantile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "quantile";
    
    return f_name;
}


const TypeSpec& Func_empiricalQuantile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
