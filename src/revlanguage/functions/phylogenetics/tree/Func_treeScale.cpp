//
//  Func_treeScale.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/6/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "Func_treeScale.h"

#include <cstddef>

#include "ModelVector.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "TreeScaleFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_treeScale::Func_treeScale( void ) : TypedFunction<TimeTree>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_treeScale* Func_treeScale::clone( void ) const {
    
    return new Func_treeScale( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_treeScale::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* scale                             = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau               = static_cast<const TimeTree&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    size_t nTips = tau->getValue().getNumberOfTips();
    RevBayesCore::RbVector<double> tipAges;
    if ( this->args[2].getVariable()->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
//        tipAges = static_cast< RevBayesCore::RbVector<RealPos>& >( this->args[2].getVariable()->getRevObject() ).getValue();
        tipAges = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode()->getValue();
        
        // sanity check
        if ( nTips != tipAges.size() )
        {
            throw RbException( "The number of tip ages does not match the number of tips" );
        }
    }
    else
    {
        double v = static_cast<RealPos&>( this->args[2].getVariable()->getRevObject() ).getValue();
        for (size_t i = 0; i < nTips; i++)
        {
            tipAges.push_back( v );
        }
    }
    
    RevBayesCore::TreeScaleFunction* f = new RevBayesCore::TreeScaleFunction( tau, scale, tipAges );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_treeScale::getArgumentRules( void ) const
{
    
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set ) {
    
        argumentRules.push_back( new ArgumentRule( "scale",       RealPos::getClassTypeSpec(), "The multiplicator by which to scale the tree,", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "tree",        TimeTree::getClassTypeSpec(), "The tree which will be re-scaled.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        std::vector<TypeSpec> tipAgeTypes;
        tipAgeTypes.push_back( RealPos::getClassTypeSpec() );
        tipAgeTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "tipAges"    , tipAgeTypes, "A vector of ages for the tips.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_treeScale::getClassType(void)
{
    
    static std::string rev_type = "Func_treeScale";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_treeScale::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_treeScale::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTreeScale";
    
    return f_name;
}


const TypeSpec& Func_treeScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
