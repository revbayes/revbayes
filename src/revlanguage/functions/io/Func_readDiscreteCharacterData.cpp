#include <cstddef>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "Func_readDiscreteCharacterData.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlBoolean.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlAminoAcidState.h"
#include "RlDnaState.h"
#include "RlNonHomologousDiscreteCharacterData.h"
#include "RlRnaState.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "RlUserInterface.h"
#include "AbstractCharacterData.h"
#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "AminoAcidState.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "NaturalNumbersState.h"
#include "NonHomologousDiscreteCharacterData.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDiscreteTaxonData.h"
#include "RlFunction.h"
#include "RnaState.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "StandardState.h" // IWYU pragma: keep

namespace RevLanguage { class RevObject; }


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readDiscreteCharacterData* Func_readDiscreteCharacterData::clone( void ) const
{
    
    return new Func_readDiscreteCharacterData( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readDiscreteCharacterData::execute( void )
{
    
    // get the information from the arguments for reading the file
    RevBayesCore::path fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    bool return_as_vector = static_cast<const RlBoolean&>( args[1].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    if ( not RevBayesCore::exists( fn ))
    {
        std::string errorStr = "";
        RevBayesCore::formatError( fn, errorStr );
        throw RbException()<<"Could not find file or path with name "<<fn;
    }
        
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<RevBayesCore::path> file_names;
    if ( RevBayesCore::is_directory( fn ))
    {
        RevBayesCore::setStringWithNamesOfFilesInDirectory( fn, file_names );
    }
    else 
    {
        file_names.push_back( fn );
    }
    
    // get the global instance of the NCL reader and clear warnings from its warnings buffer
    RevBayesCore::NclReader reader = RevBayesCore::NclReader();
    
    // the vector of matrices;
    ModelVector<AbstractHomologousDiscreteCharacterData> *m = new ModelVector<AbstractHomologousDiscreteCharacterData>();
    
    // the return value
    RevObject* retVal = NULL;
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "file_names" because some of them may not be in a format
    // that can be read.
    size_t numFilesRead = 0;
    for (auto p = file_names.begin(); p != file_names.end(); ++p)
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
        if ( my_file_type != "unknown" )
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
                ++num_matrices_read_for_this_file;
                    
                if ( (*it)->isHomologyEstablished() == true )
                {
                    if ( data_type == "DNA" )
                    {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::DnaState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::DnaState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mDNA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mDNA );
                    }
                    else if ( data_type == "RNA" )
                    {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::RnaState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::RnaState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mRNA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mRNA );
                    }
                    else if ( data_type == "Protein" )
                    {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mAA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mAA );
                    }
                    else if ( data_type == "Standard" )
                    {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mSS = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mSS );
                    }
                    else
                    {
                        --num_matrices_read_for_this_file;
                        throw RbException() << "Unknown data type \"" << data_type << "\".";
                    }
                
                }
                else
                {
                    
                    if ( data_type == "DNA" )
                    {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::DnaState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::DnaState> *>( *it );
                        NonHomologousDiscreteCharacterData<DnaState> mDNA = NonHomologousDiscreteCharacterData<DnaState>( coreM );
                        m->push_back( mDNA );
                    }
                    else if ( data_type == "RNA" )
                    {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::RnaState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::RnaState> *>( *it );
                        NonHomologousDiscreteCharacterData<RnaState> mRNA = NonHomologousDiscreteCharacterData<RnaState>( coreM );
                        m->push_back( mRNA );
                    }
                    else if ( data_type == "Protein" )
                    {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *>( *it );
                        NonHomologousDiscreteCharacterData<AminoAcidState> mAA = NonHomologousDiscreteCharacterData<AminoAcidState>( coreM );
                        m->push_back( mAA );
                    }
                    else
                    {
                        num_matrices_read_for_this_file--;
                        throw RbException() << "Unknown data type \"" << data_type << "\".";
                    }
                    
                }
                    
            }
                
        }
        else
        {
            reader.addWarning("Unknown file type");
        }
        
        if (num_matrices_read_for_this_file > 0)
        {
            ++numFilesRead;
        }
    }
    
    // print summary of results of file reading to the user
    if (RevBayesCore::is_directory( fn ))
    {
        std::stringstream o2;
        if ( numFilesRead == 0 )
        {
            o2 << "Failed to read any files from directory " << fn;
        }
        else if ( numFilesRead == 1 )
        {
            if ( m->size() == 1 )
            {
                o2 << "Successfully read one file with one character matrix from directory " << fn;
            }
            else
            {
                o2 << "Successfully read one file with " << m->size() << " character matrices from directory " << fn;
            }
        }
        else
        {
            o2 << "Successfully read " << numFilesRead << " files with " << m->size() << " character matrices from directory " << fn;
        }
        RBOUT(o2.str());
        std::set<std::string> myWarnings = reader.getWarnings();
        if ( file_names.size() - numFilesRead > 0 && myWarnings.size() > 0 )
        {
            std::stringstream o3;
            if (file_names.size() - numFilesRead == 1)
            {
                o3 << "Did not read a file for the following ";
            }
            else
            {
                    o3 << "Did not read " << file_names.size() - numFilesRead << " files for the following ";
            }
            
            if (myWarnings.size() == 1)
            {
                o3 << "reason:";
            }
            else
            {
                o3 << "reasons:";
            }
            RBOUT(o3.str());
            for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
            {
                RBOUT("* "+(*it));
            }
        }

        // set the return value
        retVal = m;
    }
    else
    {
        if (m->size() == 1)
        {
            RBOUT("Successfully read one character matrix from file '" + fn.string() + "'");
                
            // set the return value
            if ( return_as_vector == false )
            {
                retVal = new AbstractHomologousDiscreteCharacterData( (*m)[0] );
                delete m;
            }
            else
            {
                retVal = m;
            }
        }
        else if (m->size() > 1)
        {
            std::stringstream o3;
            o3 << "Successfully read " << m->size() << " character matrices from file " << fn;
            RBOUT(o3.str());
            
            // set the return value
            retVal = m;
        }
        else
        {
            std::stringstream o3;
            o3 << "Error reading file " << fn;
            RBOUT(o3.str());
            RBOUT("No data matrix was read.");
            std::set<std::string> myWarnings = reader.getWarnings();
            if ( myWarnings.size() > 0 )
            {
                for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
                {
                    RBOUT("Error:   " + (*it));
                }
        
            }
        
        }
    }
    
    return new RevVariable( retVal );
}


/** Get argument rules */
const ArgumentRules& Func_readDiscreteCharacterData::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file or directory from which to read in the character data matrix.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "alwaysReturnAsVector", RlBoolean::getClassTypeSpec(), "Should we always return the character data matrix as a vector of matrices even if there is only one?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readDiscreteCharacterData::getClassType(void)
{
    
    static std::string rev_type = "Func_readDiscreteCharacterData";
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readDiscreteCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readDiscreteCharacterData::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readDiscreteCharacterData";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readDiscreteCharacterData::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readDiscreteCharacterData::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec();
    return return_typeSpec;
}

