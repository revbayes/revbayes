//
//  Func_DECRoot.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "DispersalExtinctionRootStructureFunction.h"
#include "Func_DECRoot.h"
#include "Natural.h"
#include "RateMatrix_JC.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDagNode.h"
#include <vector>
#include <cmath>

using namespace RevLanguage;

/** default constructor */
Func_DECRoot::Func_DECRoot( void ) : TypedFunction<Simplex>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_DECRoot* Func_DECRoot::clone( void ) const
{
    
    return new Func_DECRoot( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::Simplex >* Func_DECRoot::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* rf = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode<RevBayesCore::Simplex >* rs = NULL;
    if ( this->args[1].getVariable() != NULL && this->args[1].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        rs = static_cast<const Simplex&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    }
    else {
        size_t n = log2(rf->getValue().size());
        double p = 1.0 / n;
        rs = new RevBayesCore::ConstantNode<RevBayesCore::Simplex>("", new RevBayesCore::Simplex(n,p));
    }
    
//    RevBayesCore::TypedDagNode<long>* mrs = static_cast<const Natural&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::DispersalExtinctionRootStructureFunction* f = new RevBayesCore::DispersalExtinctionRootStructureFunction( rf,rs );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_DECRoot::getArgumentRules( void ) const
{
    
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "rootFreqs", ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rangeSize", Simplex::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
//        argumentRules.push_back( new ArgumentRule( "maxRangeSize", Natural::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(RbConstants::Integer::max) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_DECRoot::getClassType(void)
{
    
    static std::string rev_type = "Func_DECRoot";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_DECRoot::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_DECRoot::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDECRoot";
    
    return f_name;
}


const TypeSpec& Func_DECRoot::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
