#include <cstddef>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "ContinuousCharacterData.h"
#include "Func_readContinuousCharacterData.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlBoolean.h"
#include "RlContinuousCharacterData.h"
#include "RlString.h"
#include "RlUserInterface.h"
#include "AbstractCharacterData.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readContinuousCharacterData* Func_readContinuousCharacterData::clone( void ) const {
    
    return new Func_readContinuousCharacterData( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readContinuousCharacterData::execute( void )
{
    
    // get the information from the arguments for reading the file
    RevBayesCore::path fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    bool return_as_vector = static_cast<const RlBoolean&>( args[1].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    if ( not RevBayesCore::exists( fn ))
    {
        std::string errorStr = "";
        RevBayesCore::formatError(fn, errorStr);
        throw RbException()<<"Could not find file or path with name "<<fn;
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<RevBayesCore::path> vector_of_file_names;
    if ( RevBayesCore::is_directory( fn ))
    {
        RevBayesCore::setStringWithNamesOfFilesInDirectory( fn, vector_of_file_names);
    }
    else
    {
        vector_of_file_names.push_back( fn );
    }
    
    // get the global instance of the NCL reader and clear warnings from its warnings buffer
    RevBayesCore::NclReader reader = RevBayesCore::NclReader();
    
    // the vector of matrices;
    ModelVector<ContinuousCharacterData> *m = new ModelVector<ContinuousCharacterData>();
    
    // the return value
    RevObject* ret_val = NULL;
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.
    size_t num_files_read = 0;
    for (auto p = vector_of_file_names.begin(); p != vector_of_file_names.end(); ++p)
    {
        bool is_interleaved = false;
        std::string my_file_type = "unknown";
        std::string data_type = "unknown";
        if (reader.isNexusFile(*p) == true)
        {
            my_file_type = "nexus";
        }
        else if (reader.isPhylipFile(*p, data_type, is_interleaved) == true)
        {
            my_file_type = "phylip";
        }
        else if (reader.isFastaFile(*p, data_type) == true)
        {
            my_file_type = "fasta";
        }
        
        int num_matrices_read_for_this_file = 0;
        if (my_file_type != "unknown")
        {
            std::string suffix = "|" + data_type;
            if ( my_file_type == "phylip" )
            {
                if (is_interleaved == true)
                {
                    suffix += "|interleaved";
                }
                else
                {
                    suffix += "|noninterleaved";
                }
            }
            else if ( my_file_type == "fasta" )
            {
                suffix += "|noninterleaved";
            }
            else
            {
                suffix += "|unknown";
            }
            my_file_type += suffix;
            
            // read the content of the file now
            std::vector<RevBayesCore::AbstractCharacterData*> m_i = reader.readMatrices( *p, my_file_type );
            for (std::vector<RevBayesCore::AbstractCharacterData*>::iterator it = m_i.begin(); it != m_i.end(); it++)
            {
                
                data_type = (*it)->getDataType();
                
                // Assume success; correct below if failure
                num_matrices_read_for_this_file++;
                
                if ( data_type == "Continuous" )
                {
                    RevBayesCore::ContinuousCharacterData *coreM = static_cast<RevBayesCore::ContinuousCharacterData *>( *it );
                    
                    ContinuousCharacterData mCC = ContinuousCharacterData (coreM );
                    m->push_back( mCC );
                }
                else
                {
                    num_matrices_read_for_this_file--;
                    throw RbException("Unknown data type \"" + data_type + "\".");
                }
            }
        }
        else
        {
            reader.addWarning("Unknown file type");
        }
        
        if (num_matrices_read_for_this_file > 0)
        {
            num_files_read++;
        }
    }
    
    
    // print summary of results of file reading to the user
    if (RevBayesCore::is_directory( fn ))
    {
        std::stringstream o2;
        if ( num_files_read == 0 )
        {
            o2 << "Failed to read any files from directory " << fn;
        }
        else if ( num_files_read == 1 )
        {
            if ( m->size() == 1 )
            {
                o2 << "Successfully read one file with one character matrix from directory " << fn;
                ret_val = new ContinuousCharacterData( (*m)[0] );
                delete m;
            }
            else
            {
                o2 << "Successfully read one file with " << m->size() << " character matrices from directory " << fn;
                ret_val = m;
            }
        }
        else
        {
            o2 << "Successfully read " << num_files_read << " files with " << m->size() << " character matrices from directory " << fn;
            ret_val = m;
        }
        RBOUT(o2.str());
        std::set<std::string> my_warnings = reader.getWarnings();
        if ( vector_of_file_names.size() - num_files_read > 0 && my_warnings.size() > 0 )
        {
            std::stringstream o3;
            if (vector_of_file_names.size() - num_files_read == 1)
            {
                o3 << "Did not read a file for the following ";
            }
            else
            {
                o3 << "Did not read " << vector_of_file_names.size() - num_files_read << " files for the following ";
            }
            
            if (my_warnings.size() == 1)
            {
                o3 << "reason:";
            }
            else
            {
                o3 << "reasons:";
            }
            RBOUT(o3.str());
            for (std::set<std::string>::iterator it = my_warnings.begin(); it != my_warnings.end(); it++)
            {
                RBOUT("* "+(*it));
            }
            
        }
    }
    else
    {
        if (m->size() == 1)
        {
            RBOUT("Successfully read one character matrix from file '" + fn.string() + "'");
            
            // set the return value
            if ( return_as_vector == false )
            {
                ret_val = new ContinuousCharacterData( (*m)[0] );
                delete m;
            }
            else
            {
                ret_val = m;
            }
        }
        else if (m->size() > 1)
        {
            std::stringstream o3;
            o3 << "Successfully read " << m->size() << " character matrices from file " << fn;
            RBOUT(o3.str());
            
            // set the return value
            ret_val = m;
        }
        else
        {
            std::stringstream o3;
            o3 << "Error reading file " << fn;
            RBOUT(o3.str());
            RBOUT("No data matrix was read.");
            std::set<std::string> my_warnings = reader.getWarnings();
            if ( my_warnings.size() > 0 )
            {
                for (std::set<std::string>::iterator it = my_warnings.begin(); it != my_warnings.end(); it++)
                {
                    RBOUT("Error:   " + (*it));
                }
            }
        }
    }
    
    return new RevVariable( ret_val );
}


/** Get argument rules */
const ArgumentRules& Func_readContinuousCharacterData::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file or directory for the character data matrices.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "alwaysReturnAsVector", RlBoolean::getClassTypeSpec(), "Should we return this object as a vector even if it is just a single matrix?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readContinuousCharacterData::getClassType(void)
{
    
    static std::string rev_type = "Func_readContinuousCharacterData";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readContinuousCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readContinuousCharacterData::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readContinuousCharacterData";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readContinuousCharacterData::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readContinuousCharacterData::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<ContinuousCharacterData>::getClassTypeSpec();
    return return_typeSpec;
}

