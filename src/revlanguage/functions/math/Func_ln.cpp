//
//  Func_exp.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/7/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "LnFunction.h"
#include "Func_ln.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_ln::Func_ln( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ln* Func_ln::clone( void ) const {
    
    return new Func_ln( *this );
}


RevBayesCore::TypedFunction<double>* Func_ln::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::LnFunction* f = new RevBayesCore::LnFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_ln::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", RealPos::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_ln::getClassType(void)
{
    
    static std::string rev_type = "Func_ln";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_ln::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ln::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "ln";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_ln::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Jordan Koch" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_ln::getHelpDescription(void) const
{
    std::string description = "";
    description += "Returns the natural log of a (positive) value.";
    
    return description;
}

/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_ln::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example =  "# create a stochastic node with an exponential distribution\n";
    example += "x ~ dnExponential(1)\n";
    example += "\n";
    example += "# create a determinstic node that takes the natural log of x\n";
    example += "y := ln(x)\n";
    example += "\n";
    example += "# print the values for x and y\n";
    example += "x # x has the stochastic value of 2.940149\n";
    example += "y # y has the determined value of 1.07846\n";

    return example;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_ln::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "log" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_ln::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Natural log function";
    
    return title;
}



const TypeSpec& Func_ln::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
