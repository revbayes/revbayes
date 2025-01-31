#include "Func_readAncestralStateTreeTrace.h"

#include <cmath>
#include <cstddef>
#include <map>
#include <iostream>
#include <string>
#include <vector>

#include "Delimiter.h"
#include "NewickConverter.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlString.h"
#include "RlTraceTree.h"
#include "StringUtilities.h"
#include "TreeUtilities.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Delimiter.h"
#include "Integer.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Trace.h"
#include "TraceTree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readAncestralStateTreeTrace* Func_readAncestralStateTreeTrace::clone( void ) const
{
    
    return new Func_readAncestralStateTreeTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readAncestralStateTreeTrace::execute( void )
{
    
    // get the information from the arguments for reading the file
    RevBayesCore::path fn        = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  treetype = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string&  sep      = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    if ( not RevBayesCore::exists( fn ))
    {
        std::string errorStr = "";
        RevBayesCore::formatError(fn, errorStr);
        throw RbException(errorStr);
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<RevBayesCore::path> vectorOfFileNames;
    if (RevBayesCore::is_directory(fn))
    {
        RevBayesCore::setStringWithNamesOfFilesInDirectory( fn, vectorOfFileNames );
    }
    else
    {
        vectorOfFileNames.push_back( fn );
    }
    
    TraceTree *rv;
    if ( treetype == "clock" )
    {
        rv = readTimeTrees(vectorOfFileNames, sep);
    }
    else if ( treetype == "non-clock" )
    {
        rv = readBranchLengthTrees(vectorOfFileNames, sep);
    }
    else
    {
        throw RbException("Unknown tree type to read.");
    }
    
    long burnin = 0;

    RevObject& b = args[3].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burninFrac = static_cast<const Probability &>(b).getValue();
        burnin = long( floor( rv->getValue().size()*burninFrac ) );
    }

    rv->getValue().setBurnin(burnin);

    return new RevVariable( rv );
}


/** Get argument rules */
const ArgumentRules& Func_readAncestralStateTreeTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
		
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        std::vector<std::string> options;
        options.push_back( "clock" );
        options.push_back( "non-clock" );
        argumentRules.push_back( new OptionRule( "treetype", new RlString("clock"), options, "The type of tree." ) );

        argumentRules.push_back( new Delimiter() );

        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes     , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );


        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readAncestralStateTreeTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_readAncestralStateTreeTrace";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readAncestralStateTreeTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readAncestralStateTreeTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readAncestralStateTreeTrace";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readAncestralStateTreeTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readAncestralStateTreeTrace::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TraceTree::getClassTypeSpec();
    return return_typeSpec;
}


TraceTree* Func_readAncestralStateTreeTrace::readBranchLengthTrees(const std::vector<RevBayesCore::path> &vectorOfFileNames, const std::string &delimiter)
{
    
    
    std::vector<RevBayesCore::TraceTree> data;
    
    for (auto& fn : vectorOfFileNames)
    {
        bool hasHeaderBeenRead = false;
        
        /* Open file */
        std::ifstream inFile( fn.string() );
        
        if ( !inFile )
            throw RbException()<< "Could not open file "<<fn;
        
        /* Initialize */
        std::string commandLine;
        std::cout << "Processing file " << fn << std::endl;
        
        size_t index = 0;
        
        std::string outgroup = "";
        
        /* Command-processing loop */
        while ( inFile.good() )
        {
            
            // Read a line
            std::string line;
            RevBayesCore::safeGetline(inFile, line);
            
            // skip empty lines
            //line = StringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
            if (line.length() == 0)
            {
                continue;
            }
            
            
            // removing comments
            if (line[0] == '#')
            {
                continue;
            }
            
            // splitting every line into its columns
            std::vector<std::string> columns;
            
            // we should provide other delimiters too
            StringUtilities::stringSplit(line, delimiter, columns);
            
            // we assume a header at the first line of the file
            if (!hasHeaderBeenRead) {
                
                for (size_t j=1; j<columns.size(); j++) {
                    
                    std::string parmName = columns[j];
                    if ( parmName == "Posterior" || parmName == "Likelihood" || parmName == "Prior") {
                        continue;
                    }
                    index = j;
                    
                    RevBayesCore::TraceTree t = RevBayesCore::TraceTree( false );
                    
                    t.setParameterName(parmName);
                    t.setFileName(fn);
                    
                    data.push_back( t );
                }
                
                hasHeaderBeenRead = true;
                
                continue;
            }
            
            // adding values to the Traces
            RevBayesCore::TraceTree& t = data[0];
            
            RevBayesCore::NewickConverter c;
            RevBayesCore::Tree *tau = c.convertFromNewick( columns[index] );
			
            t.addObject( tau );
        }
    }	
    return new TraceTree( data[0] );
}


TraceTree* Func_readAncestralStateTreeTrace::readTimeTrees(const std::vector<RevBayesCore::path> &vectorOfFileNames, const std::string &delimiter) {
    
    
    std::vector<RevBayesCore::TraceTree> data;
    
    for (auto& fn: vectorOfFileNames)
    {
        bool hasHeaderBeenRead = false;
        
        /* Open file */
        std::ifstream inFile( fn.string() );
        
        if ( !inFile )
            throw RbException()<<"Could not open file "<<fn;
        
        /* Initialize */
        std::string commandLine;
        std::cout << "Processing file \"" << fn << "\"" << std::endl;
        
        size_t index = 0;
        
        /* Command-processing loop */
        while ( inFile.good() )
        {
            
            // Read a line
            std::string line;
            RevBayesCore::safeGetline(inFile, line);
            
            // skip empty lines
            //line = StringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
            if (line.length() == 0)
            {
                continue;
            }
            
            
            // removing comments
            if (line[0] == '#')
            {
                continue;
            }
            
            // splitting every line into its columns
            std::vector<std::string> columns;
            
            // we should provide other delimiters too
            StringUtilities::stringSplit(line, delimiter, columns);
            
            // we assume a header at the first line of the file
            if ( hasHeaderBeenRead == false )
            {
                
                for (size_t j=1; j<columns.size(); j++)
                {
                    
                    std::string parmName = columns[j];
                    if ( parmName == "Posterior" || parmName == "Likelihood" || parmName == "Prior")
                    {
                        continue;
                    }
                    index = j;
                    
                    RevBayesCore::TraceTree t = RevBayesCore::TraceTree( true );
                    
                    t.setParameterName(parmName);
                    t.setFileName(fn);
                    
                    data.push_back( t );
                }
                
                hasHeaderBeenRead = true;
                
                continue;
            }
            
            // adding values to the Traces
            RevBayesCore::TraceTree& t = data[0];
            
            RevBayesCore::NewickConverter c;
            RevBayesCore::Tree *blTree = c.convertFromNewick( columns[index] );
            RevBayesCore::Tree *tau = RevBayesCore::TreeUtilities::convertTree( *blTree, false );
            
            t.addObject( tau );
			
        }
    }
    
    return new TraceTree( data[0] );
}


