//
//  Func_stateCountRateModifier.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_stateCountRateModifier.h"

#include <cstddef>

#include "StateCountRateModifierFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_stateCountRateModifier::Func_stateCountRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_stateCountRateModifier* Func_stateCountRateModifier::clone( void ) const {
    
    return new Func_stateCountRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_stateCountRateModifier::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sf = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
   
    size_t nc = static_cast<const Natural&>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::StateCountRateModifierFunction* f = new RevBayesCore::StateCountRateModifierFunction(sf, nc);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_stateCountRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        
        argumentRules.push_back( new ArgumentRule( "stateFactors", ModelVector<Real>::getClassTypeSpec(), "Power penalty (bonus) for n-1 states", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numChars", Natural::getClassTypeSpec(), "Number of characters", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_stateCountRateModifier::getClassType(void)
{
    
    static std::string revType = "Func_stateCountRateModifier";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_stateCountRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_stateCountRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnStateCountRateModifier";
    
    return f_name;
}


const TypeSpec& Func_stateCountRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}