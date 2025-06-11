#include <cstddef>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "Func_readCharacterData.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "NonHomologousDiscreteCharacterData.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlAminoAcidState.h"
#include "RlBoolean.h"
#include "RlContinuousCharacterData.h"
#include "RlNonHomologousDiscreteCharacterData.h"
#include "RlDnaState.h"
#include "RlRnaState.h"
#include "RlString.h"
#include "RlUserInterface.h"
#include "AbstractCharacterData.h"
#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "AminoAcidState.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ContinuousCharacterData.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "NaturalNumbersState.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
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
#include "WorkspaceToCoreWrapperObject.h"
#include "WorkspaceVector.h"
#include "StandardState.h" // IWYU pragma: keep

using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readCharacterData* Func_readCharacterData::clone( void ) const {
    
    return new Func_readCharacterData( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readCharacterData::execute( void ) {
        
    // get the information from the arguments for reading the file
    RevBayesCore::path fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    bool returnAsVector = static_cast<const RlBoolean&>( args[1].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    if ( not RevBayesCore::exists( fn ))
    {
        std::string errorStr = "";
        RevBayesCore::formatError( fn, errorStr );
        throw RbException()<<"Could not find file or path with name "<<fn;
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<RevBayesCore::path> vectorOfFileNames;
    if ( RevBayesCore::is_directory( fn ))
    {
        RevBayesCore::setStringWithNamesOfFilesInDirectory(fn, vectorOfFileNames);
    }
    else
    {
        vectorOfFileNames.push_back( fn );
    }
    
    // get the global instance of the NCL reader and clear warnings from its warnings buffer
    RevBayesCore::NclReader reader = RevBayesCore::NclReader();
    
    // the vector of matrices;
    WorkspaceVector<RevObject>* m = new WorkspaceVector<RevObject>();
    
    // the return value
    RevObject* retVal = NULL;
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.
    size_t numFilesRead = 0;
    for (auto p = vectorOfFileNames.begin(); p != vectorOfFileNames.end(); p++)
        {
        bool isInterleaved = false;
        std::string myFileType = "unknown";
        std::string dType = "unknown";
        if (reader.isNexusFile(*p) == true)
            {
            myFileType = "nexus";
            }
        else if (reader.isPhylipFile(*p, dType, isInterleaved) == true)
            {
            myFileType = "phylip";
            }
        else if (reader.isFastaFile(*p, dType) == true)
            {
            myFileType = "fasta";
            }
        
        int numMatricesReadForThisFile = 0;
        if (myFileType != "unknown")
            {
            std::string suffix = "|" + dType;
            if ( myFileType == "phylip" )
                {
                if (isInterleaved == true)
                    {
                    suffix += "|interleaved";
                    }
                else
                    {
                    suffix += "|noninterleaved";
                    }
                }
            else if ( myFileType == "fasta" )
                {
                suffix += "|noninterleaved";
                }
            else
                {
                suffix += "|unknown";
                }
            myFileType += suffix;
            
            // read the content of the file now
            std::vector<RevBayesCore::AbstractCharacterData*> m_i = reader.readMatrices( *p, myFileType );
            for (std::vector<RevBayesCore::AbstractCharacterData*>::iterator it = m_i.begin(); it != m_i.end(); it++)
                {
                dType = (*it)->getDataType();
                
                // Assume success; correct below if failure
                numMatricesReadForThisFile++;
                
                if ( (*it)->isHomologyEstablished() == true )
                    {
                    if ( dType == "DNA" )
                        {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::DnaState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::DnaState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mDNA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mDNA );
                        }
                    else if ( dType == "RNA" )
                        {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::RnaState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::RnaState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mRNA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mRNA );
                        }
                    else if ( dType == "Protein" )
                        {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mAA = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mAA );
                        }
                    else if ( dType == "Standard" )
                        {
                        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState> *coreM = static_cast<RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState> *>( *it );
                        AbstractHomologousDiscreteCharacterData mSS = AbstractHomologousDiscreteCharacterData( coreM );
                        m->push_back( mSS );
                        }
                    else if ( dType == "Continuous" )
                        {
                        RevBayesCore::ContinuousCharacterData *coreM = static_cast<RevBayesCore::ContinuousCharacterData *>( *it );
                        ContinuousCharacterData mCC = ContinuousCharacterData (coreM );
                        m->push_back( mCC );
                        }
                    else
                        {
                        numMatricesReadForThisFile--;
                        throw RbException() << "Unknown data type \"" << dType << "\".";
                        }
                    }
                else
                    {
                    if ( dType == "DNA" )
                        {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::DnaState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::DnaState> *>( *it );
                        NonHomologousDiscreteCharacterData<DnaState> mDNA = NonHomologousDiscreteCharacterData<DnaState>( coreM );
                        m->push_back( mDNA );
                        }
                    else if ( dType == "RNA" )
                        {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::RnaState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::RnaState> *>( *it );
                        NonHomologousDiscreteCharacterData<RnaState> mRNA = NonHomologousDiscreteCharacterData<RnaState>( coreM );
                        m->push_back( mRNA );
                        }
                    else if ( dType == "Protein" )
                        {
                        RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *coreM = static_cast<RevBayesCore::NonHomologousDiscreteCharacterData<RevBayesCore::AminoAcidState> *>( *it );
                        NonHomologousDiscreteCharacterData<AminoAcidState> mAA = NonHomologousDiscreteCharacterData<AminoAcidState>( coreM );
                        m->push_back( mAA );
                        }
                    else if ( dType == "Standard" )
                        {
                        // Can't have morphological (i.e., standard) data that don't have homology established
                        }
                    else if ( dType == "Continuous" )
                        {
                        // Can't have morphological (i.e., standard) data that don't have homology established
                        }
                    else
                        {
                        numMatricesReadForThisFile--;
                        throw RbException() << "Unknown data type \"" << dType << "\".";
                        }
                    }
                }
            }
        else
            {
            reader.addWarning("Unknown file type");
            }
        
        if (numMatricesReadForThisFile > 0)
            {
            numFilesRead++;
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
                    o2 << "Successfully read one file with " << m->size() << " character matrices from directory "<<fn;
                }
            }
        else
            {
                o2 << "Successfully read " << numFilesRead << " files with " << m->size() << " character matrices from directory "<< fn;
            }
        RBOUT(o2.str());
        std::set<std::string> myWarnings = reader.getWarnings();
        if ( vectorOfFileNames.size() - numFilesRead > 0 && myWarnings.size() > 0 )
            {
            std::stringstream o3;
            if (vectorOfFileNames.size() - numFilesRead == 1)
                {
                o3 << "Did not read a file for the following ";
                }
            else
                {
                o3 << "Did not read " << vectorOfFileNames.size() - numFilesRead << " files for the following ";
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
            if ( returnAsVector == false )
                {
                retVal = (*m)[0].clone();
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
            o3 << "Successfully read " << m->size() << " character matrices from file "<<fn;
            RBOUT(o3.str());
            
            // set the return value
            retVal = m;
            }
        else
            {
            std::set<std::string> myWarnings = reader.getWarnings();
            if ( myWarnings.size() > 0 )
                {
                std::stringstream o3;
                o3 << "Error reading file " << fn;
                RBOUT(o3.str());
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
const ArgumentRules& Func_readCharacterData::getArgumentRules( void ) const {
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "File or directory names where to find the character data.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "alwaysReturnAsVector", RlBoolean::getClassTypeSpec(), "Should the value be returned as a vector even it is only a single matrix?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readCharacterData::getClassType(void) {
    
    static std::string rev_type = "Func_readCharacterData";
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readCharacterData::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readCharacterData::getFunctionName( void ) const {

    // create a name variable that is the same for all instance of this class
    std::string f_name = "readCharacterData";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readCharacterData::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readCharacterData::getReturnType( void ) const {
    
    static TypeSpec return_typeSpec = ModelVector<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec();
    return return_typeSpec;
}

