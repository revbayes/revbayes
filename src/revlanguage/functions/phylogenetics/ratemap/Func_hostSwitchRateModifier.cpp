//
//  Function_phylogeneticDistanceRateModifier.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/24/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_hostSwitchRateModifier.h"

#include "HostSwitchRateModifierFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlTree.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_hostSwitchRateModifier::Func_hostSwitchRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_hostSwitchRateModifier* Func_hostSwitchRateModifier::clone( void ) const {
    
    return new Func_hostSwitchRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_hostSwitchRateModifier::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree >* dist  = static_cast<const Tree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* scale = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::HostSwitchRateModifierFunction* f = new RevBayesCore::HostSwitchRateModifierFunction(dist, scale);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_hostSwitchRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "Host phylogeny that influences host switch rates.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "gainScale", ModelVector<RealPos>::getClassTypeSpec(), "Scales phylogenetic distance depending on gain type.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_hostSwitchRateModifier::getClassType(void)
{
    
    static std::string revType = "Func_hostSwitchRateModifier";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_hostSwitchRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_hostSwitchRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnHostSwitchRateModifier";
    
    return f_name;
}


const TypeSpec& Func_hostSwitchRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}