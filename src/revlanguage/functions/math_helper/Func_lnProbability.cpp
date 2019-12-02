//
//  Func_lnProbability.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/12/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "Func_lnProbability.h"
#include "Probability.h"
#include "Real.h"
#include "LnProbabilityFunction.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_lnProbability::Func_lnProbability( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_lnProbability* Func_lnProbability::clone( void ) const
{
    
    return new Func_lnProbability( *this );
}


RevBayesCore::TypedFunction<double>* Func_lnProbability::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<double>* y = static_cast<RevBayesCore::StochasticNode<double>* >(x);
    RevBayesCore::LnProbabilityFunction* f = new RevBayesCore::LnProbabilityFunction( y );
    
    return f;
}

/* Get argument rules */
const ArgumentRules& Func_lnProbability::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_lnProbability::getClassType(void)
{
    
    static std::string rev_type = "Func_lnProbability";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_lnProbability::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_lnProbability::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnLnProbability";
    
    return f_name;
}


const TypeSpec& Func_lnProbability::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
