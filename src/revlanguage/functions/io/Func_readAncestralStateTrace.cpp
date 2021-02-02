/**
 * @file
 * This file contains the implementation of Func_readAncestralStateTrace.
 *
 * @brief Implementation of Func_readAncestralStateTrace
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 */

#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_readAncestralStateTrace.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlString.h"
#include "RlAncestralStateTrace.h"
#include "StringUtilities.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Trace.h"
#include "TypeSpec.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readAncestralStateTrace* Func_readAncestralStateTrace::clone( void ) const
{
    
    return new Func_readAncestralStateTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readAncestralStateTrace::execute( void ) {
    
    // get the information from the arguments for reading the file
    const std::string&  fn       = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  sep      = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    RevBayesCore::RbFileManager myFileManager( fn );
    
    if ( !myFileManager.testFile() || !myFileManager.testDirectory() )
    {
        std::string errorStr = "";
        myFileManager.formatError(errorStr);
        throw RbException(errorStr);
    }
    
    if ( !myFileManager.isFile() ) 
    {
        throw RbException("readAncestralStateTrace only takes as input a single ancestral state trace file.");
    
	} 
    else 
    {
		
		RevBayesCore::RbVector<AncestralStateTrace> traceVector;
		std::vector<RevBayesCore::AncestralStateTrace> ancestral_states = readAncestralStates(myFileManager.getFullFileName(), sep);
		for ( size_t i = 0; i < ancestral_states.size(); i++ )
		{						
			traceVector.push_back( AncestralStateTrace( ancestral_states[i] ) );			
		}
		
		WorkspaceVector<AncestralStateTrace> *theVector = new WorkspaceVector<AncestralStateTrace>( traceVector );
		
		// return a vector of traces, 1 trace for each node
		return new RevVariable( theVector );
		
	}
}



/** Get argument rules */
const ArgumentRules& Func_readAncestralStateTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
		
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "The name of the file which holds the ancestral state trace.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "separator", RlString::getClassTypeSpec(), "The separater between sampled values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readAncestralStateTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_readAncestralStateTrace";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readAncestralStateTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readAncestralStateTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readAncestralStateTrace";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readAncestralStateTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readAncestralStateTrace::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = WorkspaceVector<AncestralStateTrace>::getClassTypeSpec();
    return return_typeSpec;
}


std::vector<RevBayesCore::AncestralStateTrace> Func_readAncestralStateTrace::readAncestralStates(const std::string &fileName, const std::string &delimitter)
{
    
    RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(fileName);
    
    std::vector<RevBayesCore::AncestralStateTrace> data;
	
	bool has_header_been_read = false;

	
	/* Open file */
	std::ifstream inFile( fm.getFullFileName().c_str() );
	
	if ( !inFile )
		throw RbException( "Could not open file \"" + fileName + "\"" );
	
	/* Initialize */
	std::string commandLine;
	std::cout << "Processing file \"" << fileName << "\"" << std::endl;
	
	/* Command-processing loop */
	while ( inFile.good() )
	{
		
		// Read a line
		std::string line;
        RevBayesCore::RbFileManager reader = RevBayesCore::RbFileManager();
        reader.safeGetline(inFile, line);
		
		// skip empty lines
		if (line.length() == 0)
		{
			continue;
		}
		
		// removing comments
		if (line[0] == '#') 
        {
			continue;
		}
		
        // split every line into its columns
		std::vector<std::string> columns;
		StringUtilities::stringSplit(line, delimitter, columns);
	
		// we assume a header at the first line of the file
		if (has_header_been_read == false) 
        {
			for (size_t j = 0; j < columns.size(); j++) 
            {
				// set up AncestralStateTrace objects for each node
				RevBayesCore::AncestralStateTrace t = RevBayesCore::AncestralStateTrace();
				std::string parmName = columns[j];
				t.setParameterName(parmName);
				t.setFileName(fileName);
				data.push_back( t );
			}
			has_header_been_read = true;
		} 
        else 
        {
			for (size_t j = 0; j < columns.size(); j++) 
            {
				// add values to the AncestralStateTrace objects
				std::string anc_state = columns[j];
			    data[j].addObject( anc_state );
			}
		}
	}

	return data;
}
