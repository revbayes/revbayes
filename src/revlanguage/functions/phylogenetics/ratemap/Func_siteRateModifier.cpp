//
//  Func_siteRateModifier.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/23/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_siteRateModifier.h"

#include "SiteRateModifierFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "Natural.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_siteRateModifier::Func_siteRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_siteRateModifier* Func_siteRateModifier::clone( void ) const {
    
    return new Func_siteRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_siteRateModifier::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* rate_multipliers  = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* event_classes;
    event_classes = static_cast<const ModelVector<ModelVector<Natural> > &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<std::int64_t> >* site_classes;
    site_classes = static_cast<const ModelVector<Natural> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::SiteRateModifierFunction* f = new RevBayesCore::SiteRateModifierFunction(rate_multipliers, event_classes, site_classes);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_siteRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        argumentRules.push_back( new ArgumentRule( "rateMultipliers", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "Multiplies character i by rate-multiplier r_i.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "eventClasses", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), "Assigns event classes to rate matrix events.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "siteClasses", ModelVector<Natural>::getClassTypeSpec(), "Assigns classes to the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );

        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_siteRateModifier::getClassType(void)
{
    
    static std::string revType = "Func_siteRateModifier";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_siteRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_siteRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnSiteRateModifier";
    
    return f_name;
}


const TypeSpec& Func_siteRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}
