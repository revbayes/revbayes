//
//  Func_geographicalDistance.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "Func_geographicalDistance.h"
//#include "HyperbolicSineFunction.h"
#include "ArgumentRule.h"
#include "GeographicalDistanceFunction.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_geographicalDistance::Func_geographicalDistance( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_geographicalDistance* Func_geographicalDistance::clone( void ) const
{
    
    return new Func_geographicalDistance( *this );
}


RevBayesCore::TypedFunction<double>* Func_geographicalDistance::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* lat1 = static_cast<const RealPos&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lon1 = static_cast<const RealPos&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lat2 = static_cast<const RealPos&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lon2 = static_cast<const RealPos&>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::GeographicalDistanceFunction* f = new RevBayesCore::GeographicalDistanceFunction( lat1, lon1, lat2, lon2 );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_geographicalDistance::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "lat1", Real::getClassTypeSpec(), "The latitude for the first coordinate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lon1", Real::getClassTypeSpec(), "The longitude for the first coordinate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lat2", Real::getClassTypeSpec(), "The latitude for the second coordinate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lon2", Real::getClassTypeSpec(), "The longitude for the second coordinate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_geographicalDistance::getClassType(void)
{
    
    static std::string rev_type = "Func_geographicalDistance";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_geographicalDistance::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_geographicalDistance::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGeographicalDistance";
    
    return f_name;
}


const TypeSpec& Func_geographicalDistance::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
