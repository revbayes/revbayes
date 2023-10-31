#include <cstddef>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_ls.h"
#include "OptionRule.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "RlUserInterface.h"
#include "TypeSpec.h"
#include "Workspace.h"
#include "ArgumentRules.h"
#include "Environment.h"
#include "FunctionTable.h"
#include "Procedure.h"
#include "RbBoolean.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "StringUtilities.h"

using namespace RevLanguage;

/** Default constructor */
Func_ls::Func_ls( void ) : Procedure()
{
    
}

std::string Func_ls::makeFilteredFunctionString( const std::multimap<std::string, Function*>& table, std::string filter ) const
{
    FunctionTable filteredTable;
    
    std::stringstream s;

    std::multimap<std::string, Function*>::const_iterator it;
    
    std::stringstream t1, t2;
    for (it = table.begin(); it != table.end(); it++)
    {
        std::string name = it->first;
        
        bool keep = matchToken(name, filter);
        
        if (keep) {
            
            t1.str("");
            t1 << (*it).first << " = ";
            t2.str("");
            (*it).second->printValue( t2, true );
            s << t1.str() << StringUtilities::oneLiner( t2.str(), 75 - t1.str().size() ) << std::endl;
        }
    }
    
    return s.str();
}

bool Func_ls:: matchToken( std::string name, std::string filter ) const
{

    std::string prefix = name.substr(0,2);
    
    bool match = false;
    
    if (filter == "all") {
        match = true;
    }
    else if (filter == "function" && prefix == "fn") {
        match = true;
    }
    else if (filter == "distribution" && prefix == "dn") {
        match = true;
    }
    else if (filter == "move" && prefix == "mv") {
        match = true;
    }
    else if (filter == "monitor" && prefix == "mn") {
        match = true;
    }
    else if (filter == "analysis" && matchAnalysisToken( name )) {
        match = true;
    }
    else if (filter == "math" && matchMathToken( name )) {
        match = true;
    }
    else if (filter == "file" && matchFileToken( name )) {
        match = true;
    }
    else if (filter == "summary" && matchSummaryToken( name )) {
        match = true;
    }
    else if (filter == "utility" && matchUtilityToken( name )) {
        match = true;
    }
    else if (filter == "uncategorized" && !isTokenCategorized( name ))
    {
        match = true;
    }
    
    return match;
}

bool Func_ls::matchAnalysisToken( std::string name ) const
{
    bool match = false;
    std::set<std::string> valid_tokens;
    valid_tokens.insert("mcmc");
    valid_tokens.insert("mcmcmc");
    valid_tokens.insert("model");
    valid_tokens.insert("HillClimber");
    valid_tokens.insert("CorrespondenceAnalysis");
    valid_tokens.insert("BootstrapAnalysis");
    
    if (valid_tokens.find( name ) != valid_tokens.end()) {
        match = true;
    }
    
    return match;
}

bool Func_ls::matchMathToken( std::string name ) const
{
    bool match = false;
    std::set<std::string> valid_tokens;
    valid_tokens.insert("abs");
    valid_tokens.insert("append");
    valid_tokens.insert("ceil");
    valid_tokens.insert("choose");
    valid_tokens.insert("exp");
    valid_tokens.insert("floor");
    valid_tokens.insert("gamma");
    valid_tokens.insert("ifelse");
    valid_tokens.insert("ln");
    valid_tokens.insert("log");
    valid_tokens.insert("logistic");
    valid_tokens.insert("matrix");
    valid_tokens.insert("max");
    valid_tokens.insert("mean");
    valid_tokens.insert("median");
    valid_tokens.insert("min");
    valid_tokens.insert("normalize");
    valid_tokens.insert("power");
    valid_tokens.insert("range");
    valid_tokens.insert("rep");
    valid_tokens.insert("replicate");
    valid_tokens.insert("seq");
    valid_tokens.insert("simplex");
    valid_tokens.insert("sinh");
    valid_tokens.insert("sort");
    valid_tokens.insert("sqrt");
    valid_tokens.insert("stdev");
    valid_tokens.insert("sum");
    valid_tokens.insert("tanh");
    valid_tokens.insert("trunc");
    valid_tokens.insert("v");
    valid_tokens.insert("var");
    
    if (valid_tokens.find( name ) != valid_tokens.end()) {
        match = true;
    }
    
    return match;
}

bool Func_ls::matchFileToken( std::string name ) const
{
    bool match = false;
    std::set<std::string> valid_tokens;
    valid_tokens.insert("readAncestralStateTrace");
    valid_tokens.insert("readAncestralStateTreeTrace");
    valid_tokens.insert("readAtlas");
    valid_tokens.insert("readBranchLengthTrees");
    valid_tokens.insert("readCharacterData");
    valid_tokens.insert("readCharacterDataDelimited");
    valid_tokens.insert("readContinuousCharacterData");
    valid_tokens.insert("readDataDelimitedFile");
    valid_tokens.insert("readDiscreteCharacterData");
    valid_tokens.insert("readDistanceMatrix");
    valid_tokens.insert("readMatrix");
    valid_tokens.insert("readPoMoCountFile");
    valid_tokens.insert("readRelativeNodeAgeConstraints");
    valid_tokens.insert("readRelativeNodeAgeWeightedConstraints");
    valid_tokens.insert("readStochasticVariableTrace");
    valid_tokens.insert("readTable");
    valid_tokens.insert("readTaxonData");
    valid_tokens.insert("readTrace");
    valid_tokens.insert("readTreeTrace");
    valid_tokens.insert("readTrees");
    valid_tokens.insert("write");
    valid_tokens.insert("writeCharacterDataDelimited");
    valid_tokens.insert("writeFasta");
    valid_tokens.insert("writeNexus");
    
    if (valid_tokens.find( name ) != valid_tokens.end()) {
        match = true;
    }
    
    return match;
}

bool Func_ls::isTokenCategorized( std::string name ) const
{
    std::vector<std::string> names;
    names.push_back("utility");
    names.push_back("math");
    names.push_back("distribution");
    names.push_back("function");
    names.push_back("move");
    names.push_back("monitor");
    names.push_back("analysis");
    names.push_back("file");
    names.push_back("summary");
    
    for (size_t i = 0; i < names.size(); i++) {
        if ( matchToken(name, names[i]) ) {
            return true;
        }
    }
    
    return false;
}

bool Func_ls::matchSummaryToken( std::string name ) const
{
    bool match = false;
    std::set<std::string> valid_tokens;
    valid_tokens.insert("treeTrace");
    valid_tokens.insert("mapTree");
    valid_tokens.insert("mccTree");
    valid_tokens.insert("conTree");
    valid_tokens.insert("sumt");
    valid_tokens.insert("annotateTree");
    valid_tokens.insert("branchScoreDistance");
    
    if (valid_tokens.find( name ) != valid_tokens.end()) {
        match = true;
    }
    
    return match;
}

bool Func_ls::matchUtilityToken( std::string name ) const
{
    bool match = false;
    std::set<std::string> valid_tokens;
    
    valid_tokens.insert("clear");
    valid_tokens.insert("concatenate");
    valid_tokens.insert("exists");
    valid_tokens.insert("getOption");
    valid_tokens.insert("getwd");
    valid_tokens.insert("help");
    valid_tokens.insert("license");
    valid_tokens.insert("listFiles");
    valid_tokens.insert("listOptions");
    valid_tokens.insert("ls");
    valid_tokens.insert("module");
    valid_tokens.insert("print");
    valid_tokens.insert("printSeed");
    valid_tokens.insert("q");
    valid_tokens.insert("quit");
    valid_tokens.insert("seed");
    valid_tokens.insert("setOption");
    valid_tokens.insert("setwd");
    valid_tokens.insert("str");
    valid_tokens.insert("structure");
    valid_tokens.insert("system");
    valid_tokens.insert("time");
    valid_tokens.insert("type");
    
    if (valid_tokens.find( name ) != valid_tokens.end()) {
        match = true;
    }
    
    return match;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ls* Func_ls::clone( void ) const
{
    
    return new Func_ls( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_ls::execute( void )
{
    
    bool printAll = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();
    std::string filter = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    
    // print workspace variables
    if (filter == "all") {
        printVariables(printAll, filter);
    }
   
    // print filtered functions
    if (filter != "datatype") {
        if (filter == "all") {
            
            std::vector<std::string> names;
            names.push_back("utility");
            names.push_back("math");
            names.push_back("distribution");
            names.push_back("function");
            names.push_back("move");
            names.push_back("monitor");
            names.push_back("analysis");
            names.push_back("file");
            names.push_back("summary");
            names.push_back("uncategorized");
            
            for (size_t i = 0; i < names.size(); i++) {
                printFunctions(printAll, names[i]);
            }

        } else {
             printFunctions(printAll, filter);
        }
       
    }
    
    // print filtered types
    if (filter == "all" || filter == "datatype") {
        printTypes(printAll,filter);
    }
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_ls::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "all", RlBoolean::getClassTypeSpec(), "Should we print all variables and functions including provided ones by RevBayes?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        
        std::vector<std::string> options;
        options.push_back( "all" );
        options.push_back( "utility" );
        options.push_back( "math" );
        options.push_back( "distribution" );
        options.push_back( "function" );
        options.push_back( "datatype" );
        options.push_back( "move" );
        options.push_back( "monitor" );
        options.push_back( "analysis" );
        options.push_back( "file" );
        options.push_back( "summary" );
        options.push_back( "uncategorized" );
        
        argumentRules.push_back( new OptionRule( "filter", new RlString("all"), options, "List objects for the provided type." ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_ls::getClassType(void)
{
    
    static std::string rev_type = "Func_ls";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_ls::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ls::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "ls";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_ls::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_ls::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}


void Func_ls::printTypes( bool printAll, std::string filter ) const
{
    std::stringstream s;

    
    TypeTable types = TypeTable();
    
    if ( printAll )
    {
        const TypeTable& globalTypes = Workspace::globalWorkspace().getTypeTable();
        types.insert( globalTypes.begin(), globalTypes.end() );
    }
    
    const TypeTable& userTypes = Workspace::userWorkspace().getTypeTable();
    types.insert( userTypes.begin(), userTypes.end() );
    
    if ( !types.empty() )
    {
        RBOUT( "\n" );
        RBOUT( "Datatypes:" );
        RBOUT( "===========" );
        
        s.str("");
        TypeTable::const_iterator it;
        std::ostringstream t;
        
        for ( it = types.begin(); it != types.end(); ++it )
        {
            t.str("");
            t << (*it).first << " = ";
            t << (*it).second->getTypeSpec();
            s << StringUtilities::oneLiner( t.str(), 75 ) << std::endl;
        }
        
        s << std::endl;
        RBOUT( s.str() );
        
    }
    
    return;
}


void Func_ls::printVariables( bool printAll, std::string filter ) const
{
    
    std::stringstream s;

    VariableTable variables = VariableTable();

    if ( printAll )
    {
        const VariableTable& globalVars = Workspace::globalWorkspace().getVariableTable();
        variables.insert( globalVars.begin(), globalVars.end() );
    }

    const VariableTable& userVars = env->getVariableTable();
    variables.insert( userVars.begin(), userVars.end() );

    if ( !variables.empty() )
    {
        RBOUT( "\n" );
        RBOUT( "Variables:" );
        RBOUT( "===============" );
        
        s.str("");
        VariableTable::const_iterator it;
        std::ostringstream t1, t2;
        for ( it = variables.begin(); it != variables.end(); ++it )
        {
            t1.str("");
            t1 << (*it).first << " = ";
            t2.str("");
            (*it).second->printValue( t2, true );
            s << t1.str() << StringUtilities::oneLiner( t2.str(), 75 - t1.str().size() ) << std::endl;
        }
            
        s << std::endl;
        RBOUT( s.str() );
    }
    
    return;
}


void Func_ls::printFunctions( bool printAll, std::string filter ) const
{
    
    // we want to filter thie table prior to printing
    std::multimap<std::string, Function*> functions = FunctionTable();
    if ( printAll )
    {
        const FunctionTable& globalFuncs = Workspace::globalWorkspace().getFunctionTable();
        functions.insert( globalFuncs.begin(), globalFuncs.end() );
    }
    const FunctionTable& userFuncs = env->getFunctionTable();
    functions.insert( userFuncs.begin(), userFuncs.end() );

    // generate string of filtered functions
    std::string fn_str = makeFilteredFunctionString(functions, filter);
    
    // print filter table
    std::string title = filter;
    title[0] = toupper(title[0]);

    if ( fn_str != "" )
    {
        RBOUT( "\n" );
        RBOUT( title + ":" );
        RBOUT( "===============" );
        RBOUT( fn_str );
        
    }
    
    return;
}
