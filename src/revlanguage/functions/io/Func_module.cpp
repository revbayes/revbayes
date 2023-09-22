/**
 * @file
 * This file contains the implementation of Func_module, which is
 * the function used to read commands (module) from a file.
 *
 * @brief Implementation of Func_module
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-04 18:03:37 +0200 (Fri, 04 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface RbFunction
 * @package functions
 * @since Version 1.0, 2009-09-03
 *
 * $Id: Func_module.cpp 1485 2012-05-04 16:03:37Z hoehna $
 */

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_module.h"
#include "Module.h"
#include "ModuleSystem.h"
#include "Parser.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"
#include "ArgumentRules.h"
#include "Environment.h"
#include "Procedure.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_module::Func_module( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_module* Func_module::clone( void ) const
{
    
    return new Func_module( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_module::execute( void )
{
    
    /* Get the module */
    std::string moduleName = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    const Module& mod = ModuleSystem::getModuleSystem().getModule( moduleName );
    
    Environment *execEnv = env;
    
    if ( args[1].getVariable()->getRevObject() != RevNullObject::getInstance() )
    {
        std::string ns = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    
         execEnv = env->getChildEnvironment( ns );
    }
    
//    WorkspaceVector<RevObject> *moduleArgs = new WorkspaceVector<RevObject>();
    for (size_t i = 2; i < args.size(); ++i)
    {
//        moduleArgs->push_back( args[i].getVariable()->getRevObject() );
        if ( args[i].getLabel() != "" )
        {
            if ( !execEnv->existsVariable( args[i].getLabel() ) )
            {
                execEnv->addVariable(args[i].getLabel(), args[i].getVariable() );
            }
            
        }
        else
        {
            std::cout << "Empty ellipsis argument label.\n";
        }
    }
//    execEnv->addVariable("args", moduleArgs);
//    if ( !execEnv->existsVariable("namespace") )
//    {
//        execEnv->addVariable("namespace", new RlString(ns) );
//    }
    
    /* Initialize */
    const std::vector<std::string>& commandLines = mod.getCommandLines();
    std::string command = "";
    int lineNumber = 0;
    int result = 0;     // result from processing of last command
    RBOUT("Processing module \"" + moduleName + "\"");
    
    /* Command-processing loop */
    for ( std::vector<std::string>::const_iterator it = commandLines.begin(); it != commandLines.end(); ++it)
    {
        
        // Get a line
        const std::string& line = *it;
        lineNumber++;
        
        // If previous result was 1 (append to command), we do this
        if ( result == 1 )
            command += line;
        else
            command = line;
        
        // Process the line and record result
        result = Parser::getParser().processCommand( command, execEnv );
        if ( result == 2 ) {
            std::ostringstream msg;
            msg << "Problem processing line " << lineNumber << " in module \"" << moduleName << "\"";
            throw RbException( msg.str() );
        }
    }
    
    /* Return control */
    RBOUT("Processing of module \"" + moduleName + "\" completed");
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_module::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "Relative or absolute name of module file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "namespace", RlString::getClassTypeSpec(), "Namespace used to rescue variables from overwriting.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new Ellipsis( "Additinal variables passed into the module.", RevObject::getClassTypeSpec() ) );
        
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_module::getClassType(void)
{
    
    static std::string rev_type = "Func_module";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_module::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_module::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "module";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_module::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_module::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

