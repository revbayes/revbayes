//
//  Func_readRegionalFeatures.cpp
//  rb
//
//  Created by Fábio on 5/10/22.
//

#include <algorithm>
#include <ostream>
#include <set>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "DelimitedDataReader.h"
#include "Func_readRegionalFeatures.h"
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
Func_readRegionalFeatures* Func_readRegionalFeatures::clone( void ) const
{
    
    return new Func_readRegionalFeatures( *this );
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
RevPtr<RevVariable> Func_readRegionalFeatures::execute( void )
{

//    const ModelVector<ModelVector<RlString> >& matrix = static_cast<const ModelVector<ModelVector<RlString> >&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string& filename = static_cast<RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string& delimiter = static_cast<RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const bool& header = static_cast<RlBoolean&>( args[2].getVariable()->getRevObject() ).getValue();
    size_t header_offset = 0;
    if (header) header_offset = 1;
    
    std::cout << filename << "\n";
    
    RevBayesCore::DelimitedDataReader* rdr = new RevBayesCore::DelimitedDataReader(filename, delimiter, 0);
    const std::vector<std::vector<std::string> >&data = rdr->getChars();

    std::vector<unsigned> timeIndex;
    std::vector<unsigned> featureIndex;
    std::vector<std::string> featurePath;
    std::vector<std::string> featureRelationship;
    std::vector<std::string> featureType;

    size_t nRow = data.size();
    size_t nCol = data[0].size();

    for (int colIdx=0; colIdx < nCol; colIdx++) {
        for (int rowIdx=0; rowIdx < nRow; rowIdx++) {
            if (rowIdx > 0) {
                std::string s = data[rowIdx][colIdx];
                if (data[0][colIdx] == "time_index") {
                    int x = std::stoi(s);
                    timeIndex.push_back(x);
                }
                else if (data[0][colIdx] == "feature_index") {
                    int x = std::stoi(s);
                    featureIndex.push_back(x);
                }
                else if (data[0][colIdx] == "feature_path") {
                    featurePath.push_back( s );
                }
                else if (data[0][colIdx] == "feature_relationship") {
                    featureRelationship.push_back( s );
                }
                else if (data[0][colIdx] == "feature_type") {
                    featureType.push_back( s );
                }
            }
        }
    }

    // simple input file parsing check
    if (timeIndex.size() == 0) {
        throw RbException("Func_readRegionalFeatures: input file contained no time_index entries");
    }
    if (featureIndex.size() == 0) {
        throw RbException("Func_readRegionalFeatures: input file contained no feature_index entries");
    }
    if (featurePath.size() == 0) {
        throw RbException("Func_readRegionalFeatures: input file contained no feature_path entries");
    }
    if (featureRelationship.size() == 0) {
        throw RbException("Func_readRegionalFeatures: input file contained no feature_relationship entries");
    }
    if (featureType.size() == 0) {
        throw RbException("Func_readRegionalFeatures: input file contained no feature_type entries");
    }
    
    // index by time_index, feature_index, region_index, (region_index)
    std::map<size_t, std::map<size_t, std::vector<long> > > within_categorical;
    std::map<size_t, std::map<size_t, std::vector<double> > > within_quantitative;
    std::map<size_t, std::map<size_t, std::vector<std::vector<long> > > > between_categorical;
    std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > between_quantitative;

    // error checking
    std::set<unsigned> uniqueTimeIndex;
    std::map<std::string, std::map<std::string, std::map<size_t, std::set<size_t> > > > uniqueFeatureIndex;
    uniqueFeatureIndex["within"]["categorical"] = std::map<size_t, std::set<size_t> >();
    uniqueFeatureIndex["within"]["quantitative"]  = std::map<size_t, std::set<size_t> >();
    uniqueFeatureIndex["between"]["categorical"]  = std::map<size_t, std::set<size_t> >();
    uniqueFeatureIndex["between"]["quantitative"] = std::map<size_t, std::set<size_t> >();
    
    for (size_t i = 0; i < timeIndex.size(); i++) {
        size_t time_index = timeIndex[i];
        size_t feature_index = featureIndex[i];
        std::string feature_relationship = featureRelationship[i];
        std::string feature_type = featureType[i];
        std::string feature_path = featurePath[i];
        
        // add unique entries for error checking
        uniqueTimeIndex.insert((unsigned)time_index);

        if (uniqueFeatureIndex[feature_relationship][feature_type].find(time_index) == uniqueFeatureIndex[feature_relationship][feature_type].end()) {
            
            uniqueFeatureIndex[feature_relationship][feature_type][time_index] = std::set<size_t>();
        }
        
        uniqueFeatureIndex[feature_relationship][feature_type][time_index].insert(feature_index);
    }
    
    // check same features are present for each timeslice for each process type
    std::set<std::string> s_rel;
    s_rel.insert("within");
    s_rel.insert("between");
    std::set<std::string> s_type;
    s_type.insert("categorical");
    s_type.insert("quantitative");
    
    // within/between
    for (auto it = s_rel.begin(); it != s_rel.end(); it++) {
        // categorical/quantitative
        for (auto jt = s_type.begin(); jt != s_type.end(); jt++) {
            // time_index
            auto tmp_unique = uniqueFeatureIndex[*it][*jt];
            std::set<size_t> s1;
            std::set<size_t> s2;
            for (auto kt = uniqueTimeIndex.begin(); kt != uniqueTimeIndex.end(); kt++) {
                if (tmp_unique.find(*kt) == tmp_unique.end())
                    throw RbException() << "Missing entry in readRegionalFeatures: relationship=" << *it << " type=" << *jt << " time_index=" << *kt << "\n";
            }
            
            for (auto kt = tmp_unique.begin(); kt != tmp_unique.end(); kt++) {
                s2 = s1;
                s1 = kt->second;
                
                if (kt != tmp_unique.begin()) {
                    /*
                    std::stringstream ss;
                    ss << "r=" << *it << " t=" << *jt << " t=" << kt->first << "\n  s1=[";
                    for (std::set<size_t>::iterator ii = s1.begin(); ii != s1.end(); ii++) {
                        if (ii != s1.begin()) { ss << ","; }
                        ss << *ii;
                    }
                    ss << "]\n  s2=[";
                    for (std::set<size_t>::iterator ii = s2.begin(); ii != s2.end(); ii++) {
                        if (ii != s2.begin()) { ss << ","; }
                        ss << *ii;
                    }
                    ss << "]\n";
                    std::cout << ss.str() << "\n";
                    */
                    
                    // does set of feature_index match between time slices?
                    if (s1 != s2)
                        throw RbException() << "Missing entry in readRegionalFeatures: relationship=" << *it << " type=" << *jt << " time_index=" << kt->first << "\n";
                }
            }
        }
    }
    
    // populate
    for (size_t i = 0; i < timeIndex.size(); i++) {
        size_t time_index = timeIndex[i];
        size_t feature_index = featureIndex[i];
        std::string feature_relationship = featureRelationship[i];
        std::string feature_type = featureType[i];
        std::string feature_path = featurePath[i];
        
        RevBayesCore::DelimitedDataReader* row_rdr = new RevBayesCore::DelimitedDataReader(feature_path, delimiter, header_offset);
        std::vector<std::vector<std::string> > row_dat = row_rdr->getChars();
        
        if (feature_relationship == "within" && feature_type == "categorical") {
            for (size_t k = 0; k < row_dat[0].size(); k++) {
                long val = std::stoi( row_dat[0][k] );
                within_categorical[time_index][feature_index].push_back(val);
            }
        } else if (feature_relationship == "within" && feature_type == "quantitative") {
            for (size_t k = 0; k < row_dat[0].size(); k++) {
                double val = std::stod(row_dat[0][k] );
                within_quantitative[time_index][feature_index].push_back(val);
            }
        } else if (feature_relationship == "between" && feature_type == "categorical") {
            for (size_t j = 0; j < row_dat.size(); j++) {
                between_categorical[time_index][feature_index].push_back( std::vector<long>() );
                for (size_t k = 0; k < row_dat[0].size(); k++) {
                    long val = std::stoi( row_dat[j][k] );
                    between_categorical[time_index][feature_index][j].push_back(val);
                }
            }
        } else if (feature_relationship == "between" && feature_type == "quantitative") {
            for (size_t j = 0; j < row_dat.size(); j++) {
                between_quantitative[time_index][feature_index].push_back( std::vector<double>() );
                for (size_t k = 0; k < row_dat[0].size(); k++) {
                    double val = std::stod( row_dat[j][k] );
                    between_quantitative[time_index][feature_index][j].push_back(val);
                }
            }
        }
    }
    
    // create backend core function object
    RevBayesCore::RegionalFeatures* rf = new RevBayesCore::RegionalFeatures(within_categorical, within_quantitative, between_categorical, between_quantitative);
    
    return new RevVariable( new RlRegionalFeatures(rf) );
}


/** Get argument rules */
const ArgumentRules& Func_readRegionalFeatures::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(),
                                                  "A data table that contains the regional feature information.",
                                                  ArgumentRule::BY_VALUE,
                                                  ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "delimiter", RlString::getClassTypeSpec(), "The delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "," ) ) );
       
        argumentRules.push_back( new ArgumentRule( "header", RlBoolean::getClassTypeSpec(), "Do the summary file and the data files have headers?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
       
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readRegionalFeatures::getClassType(void)
{
    
    static std::string rev_type = "Func_readRegionalFeatures";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readRegionalFeatures::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readRegionalFeatures::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readRegionalFeatures";

    return f_name;
}


/**
 * Get the primary Rev name for this function.
 */
std::vector<std::string> Func_readRegionalFeatures::getFunctionNameAliases( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::vector<std::string> f_names;
    f_names.push_back("readTable");

    return f_names;
}


/** Get type spec */
const TypeSpec& Func_readRegionalFeatures::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readRegionalFeatures::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = WorkspaceVector<WorkspaceVector<AbstractModelObject> >::getClassTypeSpec();
    return return_typeSpec;
}




