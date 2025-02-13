#include "Func_readFossilCountsFile.h"

#include <cstdlib>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "DelimitedDataReader.h"
#include "Delimiter.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RlTaxon.h"
#include "Real.h"
#include "RealPos.h"
#include "Natural.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealPos.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "WorkspaceVector.h"
#include "AbstractModelObject.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "Integer.h"
#include "MatrixReal.h"
#include "RbBoolean.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readFossilCountsFile* Func_readFossilCountsFile::clone( void ) const
{
    
    return new Func_readFossilCountsFile( *this );
}


std::string Func_readFossilCountsFile::bitToState(const std::string &s)
{
    
    std::stringstream ss;
    char* ptr;
    long parsed = strtol(s.c_str(), &ptr, 2);
    
    if (parsed > RbConstants::Integer::max)
    {
        throw RbException("ERROR: readTSVBitsetData token " + s + " too large to store as NaturalNumber");
    }
    
    if (s.find_first_not_of("01") != std::string::npos)
    {
        throw RbException("ERROR: readTSVBitsetData token " + s + " contains non-binary characters");
    }
    
    ss << parsed;
    
    return ss.str();
}


/** Execute function */
RevPtr<RevVariable> Func_readFossilCountsFile::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string& fn                        = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    bool header                                  = static_cast<const RlBoolean&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string& del                       = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    bool rownames                                = static_cast<const RlBoolean&>( args[3].getVariable()->getRevObject() ).getValue();
    const ModelVector<Taxon>& taxa               = static_cast<const ModelVector<Taxon>&>( args[4].getVariable()->getRevObject() ).getValue();

    // get data from file
    RevBayesCore::DelimitedDataReader* tsv_data = new RevBayesCore::DelimitedDataReader(fn, del, header);
    const std::vector<std::vector<std::string> >&data = tsv_data->getChars();

    // vector to extract first column (i.e. taxa names)
    std::vector<std::string> first_elems;

    // matrix to extract fossil counts
    Natural fossil_count_matrix[ data.size() ][ data[0].size() - 1 ];
    
    // iterate through data, and fill up first_elems and fossil_count_matrix
    for (size_t i = 0; i < data.size(); ++i)
    {
        first_elems.push_back(data[i][0] );
        for (size_t j = 1; j < data[i].size(); ++j)
        {
            fossil_count_matrix[i][j - 1] = (Natural)StringUtilities::asIntegerNumber(data[i][j]);
        }
    }

    // get the sorting indices 
    std::vector<uint32_t> sorted_indices = StringUtilities::stringSortIndices(first_elems);

    // sort first_elems alphabetically
    std::sort(first_elems.begin(), first_elems.end());

    // declare resorted fossil count matrix
    Natural sorted_fossil_count_matrix[ data.size() ][ data[0].size() - 1 ];

    // iterate through, check that first elems are the same as taxa names (if a vector), and resort fossil count matrix
    for (size_t i = 0; i < data.size(); ++i)
    {
        if (std::size(first_elems) > 1 && first_elems[i] != taxa[i].getName())
        {
            throw RbException("ERROR: first column of fossil counts matrix must correspond to names in taxa vector");
        }

        for (size_t j = 1; j < data[i].size(); ++j)
        {
            sorted_fossil_count_matrix[i][j - 1] = fossil_count_matrix[sorted_indices[i]][j - 1];
        }
    }

    // create matrix return (always a natural)
    ModelVector<ModelVector<Natural> > matrix;

    // final run through to get and return matrix
    for (size_t i = 0; i < data.size(); i++)
    {
        ModelVector<Natural> r;
        for (size_t j = 1; j < data[i].size(); j++)
        {
            r.push_back(sorted_fossil_count_matrix[i][j - 1]);
        }
        matrix.push_back(r);
    }

    // return matrix
    return new RevVariable( new ModelVector<ModelVector<Natural> >(matrix) );
}


/** Get argument rules */
const ArgumentRules& Func_readFossilCountsFile::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file",      RlString::getClassTypeSpec(), "The name of the file to read.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "header",    RlBoolean::getClassTypeSpec(), "Skip first line?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        argumentRules.push_back( new Delimiter() );
        argumentRules.push_back( new ArgumentRule( "rownames",  RlBoolean::getClassTypeSpec(), "Skip first column?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        argumentRules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa corresponding to the fossil counts.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ));
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readFossilCountsFile::getClassType(void)
{
    
    static std::string rev_type = "Func_readFossilCountsFile";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readFossilCountsFile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readFossilCountsFile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readFossilCountsFile";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readFossilCountsFile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readFossilCountsFile::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = WorkspaceVector<WorkspaceVector<AbstractModelObject> >::getClassTypeSpec();
    return return_typeSpec;
}




