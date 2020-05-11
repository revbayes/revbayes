//
//  Func_exp.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/7/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <iosfwd>
#include <vector>

#include "Func_dppMeanFromConc.h"
#include "RbStatisticsHelper.h"
#include "RealPos.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_dppMeanFromConc::Func_dppMeanFromConc( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_dppMeanFromConc* Func_dppMeanFromConc::clone( void ) const
{
    
    return new Func_dppMeanFromConc( *this );
}


RevPtr<RevVariable> Func_dppMeanFromConc::execute()
{
	
    double cp = static_cast<const RealPos &>( args[0].getVariable()->getRevObject() ).getValue();
    double ne = static_cast<const RealPos &>( args[1].getVariable()->getRevObject() ).getValue();
	
	double numCat = RevBayesCore::RbStatistics::Helper::dppExpectNumTableFromConcParam(cp, ne);
	
	RealPos* value = new RealPos( numCat );
	return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_dppMeanFromConc::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "concentration", RealPos::getClassTypeSpec(), "The concentration parameter of the DPP.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numElements" , RealPos::getClassTypeSpec(), "The number of elements of the DPP.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_dppMeanFromConc::getClassType(void)
{
    
    static std::string rev_type = "Func_dppMeanFromConc";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_dppMeanFromConc::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_dppMeanFromConc::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDppMeanFromConc";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_dppMeanFromConc::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RealPos::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_dppMeanFromConc::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
