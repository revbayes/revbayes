//
//  Func_epoch.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "Func_epoch.h"
#include "EpochRateMatrixFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_epoch::Func_epoch( void ) : TypedFunction<RateGenerator>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_epoch* Func_epoch::clone( void ) const
{
    
    return new Func_epoch( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_epoch::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* t = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* r = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    
    // sanity check
    if ( t->getValue().size() != rm->getValue().size() )
    {
        throw RbException( "The number of rate matrices and epochs do not match" );
    }
    
    RevBayesCore::EpochRateMatrixFunction* f = new RevBayesCore::EpochRateMatrixFunction( rm, t, r );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_epoch::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "Q"              , ModelVector<RateGenerator>::getClassTypeSpec(), "The per epoch rate matrices", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "times"          , ModelVector<RealPos>::getClassTypeSpec(), "The times of the epochs.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rates"          , ModelVector<RealPos>::getClassTypeSpec(), "The rate multipliers per epoch.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_epoch::getClassType(void)
{
    
    static std::string rev_type = "Func_epoch";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_epoch::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_epoch::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnEpoch";
    
    return f_name;
}


const TypeSpec& Func_epoch::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
