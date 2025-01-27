#include "Func_readTrace.h"

#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Delimiter.h"
#include "Probability.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "RlTrace.h"
#include "RlUserInterface.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Integer.h"
#include "Natural.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Trace.h"
#include "TraceNumeric.h"
#include "TypeSpec.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readTrace* Func_readTrace::clone( void ) const
{
    
    return new Func_readTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readTrace::execute( void )
{

    // get the information from the arguments for reading the file
    RevBayesCore::path trace_file_name = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    // get the column delimiter
    const std::string& delimiter = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    
    
    // check that the file/path name has been correctly specified
    if ( not RevBayesCore::exists( trace_file_name ))
    {
        std::string errorStr = "";
        RevBayesCore::formatError( trace_file_name, errorStr );
        throw RbException(errorStr);
    }
        
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<RevBayesCore::path> vectorOfFileNames;
    if ( RevBayesCore::is_directory( trace_file_name) )
    {
        RevBayesCore::setStringWithNamesOfFilesInDirectory( trace_file_name, vectorOfFileNames );
    }
    else
    {
        vectorOfFileNames.push_back( trace_file_name );
    }

        
    std::vector<RevBayesCore::TraceNumeric> data;
        
    long thinning = static_cast<const Natural&>( args[3].getVariable()->getRevObject() ).getValue();

    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.

    for (auto& filename: vectorOfFileNames)
    {
        bool hasHeaderBeenRead = false;
        
        /* Open file */
        std::ifstream inFile( filename.string() );
        
        if ( !inFile )
            throw RbException()<<"Could not open file "<<filename;
            
        /* Initialize */
        std::string commandLine;
        RBOUT("Processing file \"" + filename.string() + "\"");
        size_t n_samples = 0;
            
        /* Command-processing loop */
        while ( inFile.good() )
        {
                
            // Read a line
            std::string line;
            RevBayesCore::safeGetline(inFile, line);
                
            // skip empty lines
            //line = stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
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
            StringUtilities::stringSplit(line, delimiter, columns);
                
            // we assume a header at the first line of the file
            if (!hasHeaderBeenRead)
            {
                    
                for (size_t j=0; j<columns.size(); j++)
                {
                    RevBayesCore::TraceNumeric t;
                        
                    std::string parmName = columns[j];
                    t.setParameterName(parmName);
                    t.setFileName( filename );
                        
                    data.push_back( t );
                }
                    
                hasHeaderBeenRead = true;
                    
                continue;
            }
            
            
            // increase our sample counter
            ++n_samples;
            
            // we need to check if we skip this sample in case of thinning.
            if ( (n_samples-1) % thinning > 0 )
            {
                continue;
            }
            
            // adding values to the Tracess
            for (size_t j=0; j<columns.size(); j++)
            {
                RevBayesCore::TraceNumeric& t = static_cast<RevBayesCore::TraceNumeric&>( data[j] );
                std::string tmp = columns[j];
                double d = atof( tmp.c_str() );
                t.addObject(d);
            }
        }
    }
    
    RevObject& b = args[2].getVariable()->getRevObject();

    WorkspaceVector<Trace> *rv = new WorkspaceVector<Trace>();
    for (std::vector<RevBayesCore::TraceNumeric>::iterator it = data.begin(); it != data.end(); ++it)
    {
        int burnin = 0;

        if ( b.isType( Integer::getClassTypeSpec() ) )
        {
            burnin = (int)static_cast<const Integer &>(b).getValue();
        }
        else
        {
            double burninFrac = static_cast<const Probability &>(b).getValue();
            burnin = int( floor( it->size()*burninFrac ) );
        }

        it->computeStatistics();

        Trace trace( *it );
        trace.getValue().setBurnin(burnin);

        rv->push_back( trace );
    }
    
    // return the vector of traces
    return new RevVariable( rv );
}


/** Get argument rules */
const ArgumentRules& Func_readTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "Name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Delimiter() );

        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes     , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        argumentRules.push_back( new ArgumentRule( "thinning", Natural::getClassTypeSpec(), "The frequency of samples to read, i.e., we will only used every n-th sample where n is defined by this argument.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( 1l ) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_readTrace";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTrace";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readTrace::getReturnType( void ) const
{
    static TypeSpec return_typeSpec = WorkspaceVector<Trace>::getClassTypeSpec();
    return return_typeSpec;
}


