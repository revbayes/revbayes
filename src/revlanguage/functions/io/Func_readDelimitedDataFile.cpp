#include "Func_readDelimitedDataFile.h"

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
Func_readDelimitedDataFile* Func_readDelimitedDataFile::clone( void ) const
{
    
    return new Func_readDelimitedDataFile( *this );
}


std::string Func_readDelimitedDataFile::bitToState(const std::string &s)
{
    
    std::stringstream ss;
    char* ptr;
    std::int64_t parsed = strtol(s.c_str(), &ptr, 2);
    
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
RevPtr<RevVariable> Func_readDelimitedDataFile::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string& fn  = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    bool header            = static_cast<const RlBoolean&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string& del = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    bool rownames          = static_cast<const RlBoolean&>( args[3].getVariable()->getRevObject() ).getValue();
    
    // get data from file
    RevBayesCore::DelimitedDataReader* tsv_data = new RevBayesCore::DelimitedDataReader(fn, del, header);
    const std::vector<std::vector<std::string> >&data = tsv_data->getChars();

    WorkspaceVector<WorkspaceVector<AbstractModelObject> > matrix;

    enum Datatype { UNBOUNDED   = 0b0001,
                    UNCOUNTABLE = 0b0010,
                    NUMERIC     = 0b0100,
                    ALPHA       = 0b1000,

                    NATURAL     = 0x4,
                    INTEGER     = 0x5,
                    REALPOS     = 0x6,
                    REAL        = 0x7,
                    STRING      = 0x8
                  };

    int matrix_type = 0;

    for (size_t i = 0; i < data.size(); ++i)
    {
        WorkspaceVector<AbstractModelObject> row;

        // make row
        for (size_t j= rownames; j < data[i].size(); ++j)
        {
            int elemtype = 0;

            if ( StringUtilities::isNumber(data[i][j]) )
            {
                elemtype |= NUMERIC;

                // integer
                if ( StringUtilities::isIntegerNumber(data[i][j]) )
                {
                    int val = atoi(data[i][j].c_str() );

                    // negative integer
                    if ( val < 0)
                    {
                        elemtype |= UNBOUNDED;
                        row.push_back( Integer( val ) );
                    }
                    // positive integer
                    else
                    {
                        row.push_back( Natural( val ) );
                    }
                }
                // real
                else
                {
                    elemtype |= UNCOUNTABLE;
                    double val = atof(data[i][j].c_str() );

                    // negative real
                    if ( val < 0)
                    {
                        elemtype |= UNBOUNDED;
                        row.push_back( Real( val ) );
                    }
                    // positive real
                    else
                    {
                        row.push_back( RealPos( val ) );
                    }
                }
            }
            // string
            else
            {
                elemtype = STRING;
                row.push_back( RlString( data[i][j]) );
            }

            matrix_type = matrix_type | elemtype;
        }

        matrix.push_back(row);
    }

    if (matrix_type <= STRING)
    {
        if (matrix_type == REALPOS)
        {

            RevBayesCore::MatrixReal m = RevBayesCore::MatrixReal(matrix.size(),matrix[0].size());
            for (size_t i = 0; i<matrix.size(); i++)
            {

                for (size_t j = 0; j < matrix[i].size(); j++)
                {
                    if ( matrix[i][j].isType(RealPos::getClassTypeSpec()) == false )
                    {
                        RealPos *tmp = (RealPos*)(matrix[i][j].convertTo(RealPos::getClassTypeSpec()));
                        m[i][j] = (tmp)->getValue();
                        delete tmp;
                    }
                    else
                    {
                        m[i][j] = dynamic_cast<const RealPos&>(matrix[i][j]).getValue();
                    }
                }
            }

            return new RevVariable( new MatrixRealPos(m) );
        }
        else if (matrix_type == REAL)
        {

            RevBayesCore::MatrixReal m = RevBayesCore::MatrixReal(matrix.size(),matrix[0].size());
            for (size_t i = 0; i<matrix.size(); i++)
            {

                for (size_t j = 0; j < matrix[i].size(); j++)
                {
                    if ( matrix[i][j].isType(Real::getClassTypeSpec()) == false )
                    {
                        
                        Real *tmp = (Real*)(matrix[i][j].convertTo(Real::getClassTypeSpec()));
                        m[i][j] = (tmp)->getValue();
                        delete tmp;
                    }
                    else
                    {
                        m[i][j] = dynamic_cast<const Real&>(matrix[i][j]).getValue();
                    }
                }

            }

            return new RevVariable( new MatrixReal(m) );
        }
        else if (matrix_type == NATURAL)
        {
            ModelVector<ModelVector<Natural> > m;
            for (size_t i = 0; i<matrix.size(); i++)
            {
                ModelVector<Natural> r;
                for (size_t j = 0; j < matrix[i].size(); j++)
                {
                    r.push_back(matrix[i][j]);
                }
                m.push_back(r);
            }
            return new RevVariable( new ModelVector<ModelVector<Natural> >(m) );
        }
        else if (matrix_type == INTEGER)
        {
            ModelVector<ModelVector<Integer> > m;
            for (size_t i = 0; i<matrix.size(); i++)
            {
                ModelVector<Integer> r;
                for (size_t j = 0; j < matrix[i].size(); j++)
                {
                    r.push_back(matrix[i][j]);
                }
                m.push_back(r);
            }
            return new RevVariable( new ModelVector<ModelVector<Integer> >(m) );
        }
        else if (matrix_type == STRING)
        {
            ModelVector<ModelVector<RlString> > m;
            for (size_t i = 0; i<matrix.size(); i++)
            {
                ModelVector<RlString> r;
                for (size_t j = 0; j < matrix[i].size(); j++)
                {
                    r.push_back(matrix[i][j]);
                }
                m.push_back(r);
            }
            return new RevVariable( new ModelVector<ModelVector<RlString> >(m) );
        }
    }

    return new RevVariable( new WorkspaceVector<WorkspaceVector<AbstractModelObject> >(matrix) );
}


/** Get argument rules */
const ArgumentRules& Func_readDelimitedDataFile::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file",      RlString::getClassTypeSpec(), "The name of the file to read.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "header",    RlBoolean::getClassTypeSpec(), "Skip first line?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        argumentRules.push_back( new Delimiter() );
        argumentRules.push_back( new ArgumentRule( "rownames",  RlBoolean::getClassTypeSpec(), "Skip first column?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readDelimitedDataFile::getClassType(void)
{
    
    static std::string rev_type = "Func_readDelimitedDataFile";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readDelimitedDataFile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readDelimitedDataFile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readDelimitedDataFile";

    return f_name;
}


/**
 * Get the primary Rev name for this function.
 */
std::vector<std::string> Func_readDelimitedDataFile::getFunctionNameAliases( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::vector<std::string> f_names;
    f_names.push_back("readTable");
    f_names.push_back("readDataDelimitedFile");

    return f_names;
}


/** Get type spec */
const TypeSpec& Func_readDelimitedDataFile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readDelimitedDataFile::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = WorkspaceVector<WorkspaceVector<AbstractModelObject> >::getClassTypeSpec();
    return return_typeSpec;
}




