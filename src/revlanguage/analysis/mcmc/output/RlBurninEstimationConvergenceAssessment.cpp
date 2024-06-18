#include "RlBurninEstimationConvergenceAssessment.h"

#include <cstdlib>
#include <cstddef>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BurninEstimatorContinuous.h"
#include "ConvergenceDiagnosticContinuous.h"
#include "Delimiter.h"
#include "EssMax.h"
#include "EssTest.h"
#include "GelmanRubinTest.h"
#include "GewekeTest.h"
#include "HeidelbergerWelchTest.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlString.h"
#include "RlUserInterface.h"
#include "SemMin.h"
#include "StationarityTest.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "ModelVector.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "Trace.h"
#include "TraceNumeric.h"
#include "TypedDagNode.h"
#include "WorkspaceObject.h"


using namespace RevLanguage;

BurninEstimationConvergenceAssessment::BurninEstimationConvergenceAssessment() : WorkspaceObject(),
    delimiter( "\t" ),
    filenames(),
    burninMethod( "ESS" ),
    verbose( false )
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
    
    std::vector<std::string> options;
    options.push_back( "ESS" );
    options.push_back( "SEM" );
    
    ArgumentRules* burninMethodArgRules = new ArgumentRules();
    burninMethodArgRules->push_back( new OptionRule("method", options, "The burnin estimation method." ) );
    methods.addFunction( new MemberProcedure( "setBurninMethod", RlUtils::Void, burninMethodArgRules) );
    
    ArgumentRules* verboseArgRules = new ArgumentRules();
    verboseArgRules->push_back( new ArgumentRule("x", RlBoolean::getClassTypeSpec(), "Should the output be verbose?", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "verbose", RlUtils::Void, verboseArgRules) );

    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
BurninEstimationConvergenceAssessment* BurninEstimationConvergenceAssessment::clone(void) const
{
    
	return new BurninEstimationConvergenceAssessment(*this);
}



/* Map calls to member methods */
RevPtr<RevVariable> BurninEstimationConvergenceAssessment::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    RevPtr<RevVariable> retVar;
    
    if (name == "run")
    {
        found = true;
        
        bool passed = true;
        
        RBOUT("\n\t*********************************************");
        RBOUT("\tBurn-in Estimation and Convergence Assessment");
        RBOUT("\t*********************************************\n\n");
        
        RevBayesCore::BurninEstimatorContinuous *burninEst = NULL;
        
        if ( burninMethod == "ESS" )
        {
            burninEst = new RevBayesCore::EssMax();
        }
        else if ( burninMethod == "SEM" )
        {
            burninEst = new RevBayesCore::SemMin();
        }
        else
        {
            throw RbException("Unknown burnin estimation method");
        }
        
        RevBayesCore::ConvergenceDiagnosticContinuous *essTest = new RevBayesCore::EssTest();
        RevBayesCore::ConvergenceDiagnosticContinuous *gewekeTest = new RevBayesCore::GewekeTest();
        RevBayesCore::ConvergenceDiagnosticContinuous *gelmanRubinTest = new RevBayesCore::GelmanRubinTest();
        RevBayesCore::ConvergenceDiagnosticContinuous *heidelbergerTest = new RevBayesCore::HeidelbergerWelchTest();
        RevBayesCore::ConvergenceDiagnosticContinuous *stationarityTest = new RevBayesCore::StationarityTest();
        
        // read the traces
        
        // set up a vector of strings containing the name or names of the files to be read
        std::vector<RevBayesCore::path> vectorOfFileNames;
        
        for (auto& fn : filenames)
        {
            if ( is_regular_file( fn ))
            {
                vectorOfFileNames.push_back( fn );
            }
            else if ( is_directory( fn ) )
            {
                RevBayesCore::setStringWithNamesOfFilesInDirectory( fn, vectorOfFileNames );
            }
        }
        
        
        std::stringstream tmp;
        tmp << "Processing " << vectorOfFileNames.size() << ( vectorOfFileNames.size() > 1 ? " files ..." : " file ...");
        RBOUT( tmp.str() );
        
        
        
        RBOUT("\n\t-----------------------------------");
        RBOUT("\tSingle Chain Convergence Assessment");
        RBOUT("\t-----------------------------------\n\n");
        
        std::vector< std::vector<RevBayesCore::TraceNumeric> > runs(vectorOfFileNames.size(), std::vector<RevBayesCore::TraceNumeric>() );
        std::vector< size_t > burnins;
        
        for (size_t p = 0; p < vectorOfFileNames.size(); p++)
        {
            
            auto &fn = vectorOfFileNames[p];
            
            RBOUT("\tProcessing file '" + fn.string() + "'");
            
            // read in the traces from this file
            readTrace(fn, runs[p]);
            
            // add the traces to our runs
            std::vector<RevBayesCore::TraceNumeric>& data = runs[p];
            
            size_t maxBurnin = 0;
            
            // find the max burnin
            for ( size_t i = 0; i < data.size(); ++i)
            {
                size_t b = burninEst->estimateBurnin( data[i] );

                if ( maxBurnin < b )
                {
                    maxBurnin = b;
                }
            }
            
            bool failed = false;
            size_t numFailedParams = 0;
            for ( size_t i = 0; i < data.size(); ++i)
            {
                data[i].setBurnin( maxBurnin );
                data[i].computeStatistics();
                
                bool gewekeStat = gewekeTest->assessConvergence( data[i] );
                bool essStat = essTest->assessConvergence( data[i] );
//                bool gelmanStat = gelmanRubinTest->assessConvergence( data[i] );
                bool stationarityStat = stationarityTest->assessConvergence( data[i] );
                bool heidelbergerStat = heidelbergerTest->assessConvergence( data[i] );
                bool failedParam = !gewekeStat || !stationarityStat || !heidelbergerStat || !essStat;
                
                if ( failedParam == true )
                {
                    numFailedParams++;
                }
                
                failed |= failedParam;
                
                if ( verbose == true )
                {
                    RBOUT("\t\tResults for parameter '" + data[i].getParameterName() + "'\n" );
                    std::stringstream ss("");
                    ss << "\t\t\tESS = " << data[i].getESS();
                    RBOUT( ss.str() );
                    std::string p = (gewekeStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed Geweke test:\t\t\t\t" + p);
                    p = (essStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed ESS test:\t\t\t\t" + p);
//                    p = (gelmanStat ? "TRUE" : "FALSE");
//                    RBOUT("\t\t\tPassed Gelman-Rubin test:\t\t\t" + p);
                    p = (stationarityStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed Stationarity test:\t\t\t" + p);
                    p = (heidelbergerStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed Heideberger-Welch test:\t\t" + p);
                }
                
            }
            
            if ( failed )
            {
                std::stringstream ss("");
                ss << "\tConvergence assessment detected " << numFailedParams << " possible issues in file '" + fn.string() + "'.\n\n";
                RBOUT( ss.str() );
                
            }
            else
            {
                RBOUT("No failure to converge was detected in file '"+ fn.string() +"'.\n\n");
            }
            
            passed &= !failed;
        }
        
        
        // now, compare the different runs
        if ( runs.size() > 1 )
        {
            RBOUT("\n\t----------------------------------");
            RBOUT("\tMulti Chain Convergence Assessment");
            RBOUT("\t----------------------------------\n\n");
            
            std::vector<RevBayesCore::TraceNumeric>& run = runs[0];
            
            bool failed = false;
            size_t numFailedParams = 0;
            
            for (size_t j=0; j<run.size(); ++j)
            {
                
                RevBayesCore::TraceNumeric& t = run[j];
                const std::string &traceName = t.getParameterName();
                std::vector<RevBayesCore::TraceNumeric> v;
                v.push_back( t );
                
                for (size_t i=1; i<runs.size(); ++i)
                {
                    
                    size_t index = runs[i].size();
                    for (size_t k=0; k<runs[i].size(); ++k)
                    {
                        if ( runs[i][k].getParameterName() == traceName )
                        {
                            index = k;
                            break;
                        }
                    }
                    
                    if ( index == runs[i].size() )
                    {
                        throw RbException()<<"Could not find a trace for parameter '"<<traceName<<"' in file "<<runs[i][0].getFileName()<<".";
                    }
                    RevBayesCore::TraceNumeric& nextTrace = runs[i][index];
                    v.push_back( nextTrace );
            
                }
                
                
//                bool gewekeStat = gewekeTest->assessConvergence( v );
//                bool essStat = essTest->assessConvergence( v );
                bool gelmanStat = gelmanRubinTest->assessConvergence( v );
                bool stationarityStat = stationarityTest->assessConvergence( v );
//                bool heidelbergerStat = heidelbergerTest->assessConvergence( v );
                bool failedParam =  !gelmanStat || !stationarityStat;
                
                if ( failedParam == true )
                {
                    numFailedParams++;
                }
                
                failed |= failedParam;
                
                if ( verbose == true )
                {
                    RBOUT("\t\tResults for parameter '" + traceName + "'\n" );
                    std::string p = (gelmanStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed Gelman-Rubin test:\t\t\t" + p);
                    p = (stationarityStat ? "TRUE" : "FALSE");
                    RBOUT("\t\t\tPassed Stationarity test:\t\t\t" + p);
                }
            
            }
            
            
            
            if ( failed )
            {
                std::stringstream ss("");
                ss << "\tConvergence assessment detected " << numFailedParams << " possible issues.\n\n";
                RBOUT( ss.str() );
                
            }
            else
            {
                RBOUT("No failure to convergence could be detected.\n\n");
            }
            
        }
        
        RBOUT("\n");
        
        retVar = new RevVariable( new RlBoolean( passed ) );
        
    }
    else if (name == "setBurninMethod")
    {
        found = true;
        
        // get the member with give index
        burninMethod = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
        
    }
    else if (name == "verbose")
    {
        found = true;
        
        // get the member with give index
        verbose = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();
        
    }
    else
    {
        retVar = RevObject::executeMethod( name, args, found );
    }
    
    return retVar;
}


/** Get Rev type of object */
const std::string& BurninEstimationConvergenceAssessment::getClassType(void)
{
    
    static std::string rev_type = "BurninEstimationConvergenceAssessment";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& BurninEstimationConvergenceAssessment::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string BurninEstimationConvergenceAssessment::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "beca";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& BurninEstimationConvergenceAssessment::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> filenameTypes;
        filenameTypes.push_back( RlString::getClassTypeSpec() );
        filenameTypes.push_back( ModelVector<RlString>::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule("filename", filenameTypes, "The name of the file with the parameter samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new Delimiter() );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& BurninEstimationConvergenceAssessment::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void BurninEstimationConvergenceAssessment::printValue(std::ostream &o, bool user) const
{
    
    o << "BurninEstimationConvergenceAssessment";
}



void BurninEstimationConvergenceAssessment::readTrace(const RevBayesCore::path &fn, std::vector<RevBayesCore::TraceNumeric> &data)
{
    
    bool hasHeaderBeenRead = false;
    
    /* Open file */
    std::ifstream inFile( fn.string() );
    
    if ( !inFile )
        throw RbException()<<"Could not open file "<<fn;
    
    /* Initialize */
    std::string commandLine;
    
    size_t startIndex = 0;
    
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
            
            // do not add the iteration number as a trace
            if ( columns[0] == "Iteration" )
            {
                startIndex = 1;
            }
            
            for (size_t j=startIndex; j<columns.size(); j++)
            {
                RevBayesCore::TraceNumeric t;
                
                std::string parmName = columns[j];
                t.setParameterName(parmName);
                t.setFileName(fn);
                
                data.push_back( t );
            }
            
            hasHeaderBeenRead = true;
            
            continue;
        }
        
        // adding values to the Tracess
        for (size_t j=startIndex; j<columns.size(); j++)
        {
            RevBayesCore::TraceNumeric& t = data[j-startIndex];
            std::string tmp = columns[j];
            double d = atof( tmp.c_str() );
            t.addObject(d);
        }
        
    }
    
}


/** Set a member variable */
void BurninEstimationConvergenceAssessment::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "filename")
    {
        // empty the files names
        if ( var->getRevObject().getTypeSpec().isDerivedOf( RlString::getClassTypeSpec() ) )
        {
            filenames.insert( static_cast<const RlString&>( var->getRevObject() ).getValue() );
        }
        else
        {
            const std::vector<std::string> &fn = static_cast<const ModelVector<RlString> &>( var->getRevObject() ).getValue();
            for (std::vector<std::string>::const_iterator it=fn.begin(); it!=fn.end(); ++it)
            {
                filenames.insert( *it );
            }
        }
        
    }
    else if ( name == "delimiter" || name == "separator" )
    {
        delimiter = static_cast<const RlString&>( var->getRevObject() ).getValue();
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}
