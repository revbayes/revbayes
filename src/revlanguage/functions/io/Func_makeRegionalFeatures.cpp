//
//  Func_makeRegionalFeatures.cpp
//  rb
//
//  Created by Fábio on 5/10/22.
//

#include <stdlib.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "DelimitedDataReader.h"
#include "Func_makeRegionalFeatures.h"
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
#include "RlRegionalFeatures.h"
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
Func_makeRegionalFeatures* Func_makeRegionalFeatures::clone( void ) const
{
    
    return new Func_makeRegionalFeatures( *this );
}


std::string Func_makeRegionalFeatures::bitToState(const std::string &s)
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
/* This function will:
        1. store the String[][] argument as a local matrix variable;
           rows in the matrix tell us where to look for regional feature filesets
           and what kind of regional features are stored in each fileset
        2. parse the matrix values to handle these filesets
        3. construct a RevBayesCore::RegionalFeatures based on the fileset info
        4. populate the RevBayesCore::RegionalFeatures object with the relevant filesets
           that were specified in the original matrix (Step 1)
        5. return an object to the RevBayes language layer of type RevLanguage::RlRegionalFeatures
           that is a wrapper for the populated RevBayesCore::RegionalFeatures object
  */
RevPtr<RevVariable> Func_makeRegionalFeatures::execute( void )
{

    const ModelVector<ModelVector<RlString> >& matrix = static_cast<const ModelVector<ModelVector<RlString> >&>( args[0].getVariable()->getRevObject() ).getValue();
    
//    feature_name,within_or_between,categorical_or_quantitative
//    cw_feature1,within,categorical
//    qw_feature1,within,quantitative
//    cb_feature1,between,categorical
//    qb_feature1,between,quantitative

    std::vector<std::string> featureName;
    std::vector<std::string> featureRelationship;
    std::vector<std::string> featureType;
    
    size_t nRow = matrix.size();
    std::cout << nRow << "\n";
    size_t nCol = matrix[0].size();
    
//    std::vector<std::string> matrixHeaders;
    //= { "featureName", "featureRelationship", "featureType" };
    
    for (int colIdx=0; colIdx < nCol; colIdx++) {
        for (int rowIdx=0; rowIdx < nRow; rowIdx++) {
            std::string s = matrix[rowIdx][colIdx];
            std::cout << s << "\n";
            if (rowIdx > 0) {
                if (matrix[0][colIdx] == "feature_path") {
                    featureName.push_back( s );
                }
                else if (matrix[0][colIdx] == "within_or_between") {
                    featureRelationship.push_back( s );
                }
                else if (matrix[0][colIdx] == "categorical_or_quantitative") {
                    featureType.push_back( s );
                }
            }
        }
    }
    
    RevBayesCore::RegionalFeatures* rf = new RevBayesCore::RegionalFeatures(featureName, featureRelationship, featureType);
//    return new RevVariable( new RlRegionalFeatures(rf) );
    return new RevVariable( new RlRegionalFeatures(rf) );
    // return new RevVariable( new ModelVector<ModelVector<RlString> >(m) );
}


/** Get argument rules */
const ArgumentRules& Func_makeRegionalFeatures::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "feature_table",
                                                  ModelVector<ModelVector<RlString> >::getClassTypeSpec(),
                                                  "2D vector summarizing regional feature attributes.",
                                                  ArgumentRule::BY_VALUE,
                                                  ArgumentRule::ANY ) );
//        argumentRules.push_back( new ArgumentRule( "header",    RlBoolean::getClassTypeSpec(), "Skip first line?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
//        argumentRules.push_back( new ArgumentRule( "delimiter", RlString::getClassTypeSpec(), "The delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "\t" ) ) );
//        argumentRules.push_back( new ArgumentRule( "rownames",  RlBoolean::getClassTypeSpec(), "Skip first column?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_makeRegionalFeatures::getClassType(void)
{
    
    static std::string rev_type = "Func_makeRegionalFeatures";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_makeRegionalFeatures::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_makeRegionalFeatures::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "makeRegionalFeatures";

    return f_name;
}


/**
 * Get the primary Rev name for this function.
 */
std::vector<std::string> Func_makeRegionalFeatures::getFunctionNameAliases( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::vector<std::string> f_names;
    f_names.push_back("readTable");

    return f_names;
}


/** Get type spec */
const TypeSpec& Func_makeRegionalFeatures::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_makeRegionalFeatures::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = WorkspaceVector<WorkspaceVector<AbstractModelObject> >::getClassTypeSpec();
    return return_typeSpec;
}




