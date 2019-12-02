//
//  Func_mrcaIndex.cpp
//  RevLanguage
//
//  Created by Michael Landis on 8/19/14.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Clade.h"
#include "Func_mrcaIndex.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "MrcaIndexStatistic.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_mrcaIndex::Func_mrcaIndex( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_mrcaIndex* Func_mrcaIndex::clone( void ) const
{
    
    return new Func_mrcaIndex( *this );
}


RevBayesCore::TypedFunction<long>* Func_mrcaIndex::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::Clade& c = static_cast<const Clade &>( this->args[1].getVariable()->getRevObject() ).getValue();
    RevBayesCore::MrcaIndexStatistic* f = new RevBayesCore::MrcaIndexStatistic( tau, c );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_mrcaIndex::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree" , TimeTree::getClassTypeSpec(), "The tree which is used to compute the MRCA.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "clade", Clade::getClassTypeSpec()   , "The clade for which the MRCA is searched.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_mrcaIndex::getClassType(void)
{
    
    static std::string rev_type = "Func_mrcaIndex";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_mrcaIndex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_mrcaIndex::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mrcaIndex";
    
    return f_name;
}


const TypeSpec& Func_mrcaIndex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
