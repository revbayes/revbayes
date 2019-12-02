#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Ellipsis.h"
#include "Func_readTreeTrace.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "NewickConverter.h"
#include "OptionRule.h"
#include "Probability.h"
#include "ProgressBar.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlBranchLengthTree.h"
#include "RlClade.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "RlUserInterface.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "TraceTree.h"
#include "TreeUtilities.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readTreeTrace* Func_readTreeTrace::clone( void ) const
{
    
    return new Func_readTreeTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readTreeTrace::execute( void )
{
    
    size_t arg_index_files     = 0;
    size_t arg_index_tree_type = 1;
    size_t arg_index_outgroup  = 2;
    size_t arg_index_separator = 3;
    size_t arg_index_burnin    = 4;
    size_t arg_index_thinning  = 5;
    size_t arg_index_nexus  = 6;

    // get the information from the arguments for reading the file
    const std::string&  treetype = static_cast<const RlString&>( args[arg_index_tree_type].getVariable()->getRevObject() ).getValue();
    const std::string&  sep      = static_cast<const RlString&>( args[arg_index_separator].getVariable()->getRevObject() ).getValue();
    long                thin     = static_cast<const Natural&>( args[arg_index_thinning].getVariable()->getRevObject() ).getValue();
    bool nexus = static_cast<RlBoolean&>(args[arg_index_nexus].getVariable()->getRevObject()).getValue();

    std::vector<std::string> vectorOfFileNames;
    
    if ( args[0].getVariable()->getRevObject().isType( ModelVector<RlString>::getClassTypeSpec() ) )
    {
        const ModelVector<RlString>& fn = static_cast<const ModelVector<RlString> &>( args[arg_index_files].getVariable()->getRevObject() ).getValue();
        std::vector<std::string> tmp;
        for (size_t i = 0; i < fn.size(); i++)
        {
            RevBayesCore::RbFileManager myFileManager( fn[i] );
            
            if ( !myFileManager.testFile() && !(myFileManager.testDirectory() && myFileManager.getFileName() == "" ) )
            {
                std::string errorStr = "Could not find filename: " + myFileManager.getFileName() + "\n";
                myFileManager.formatError(errorStr);
                throw RbException(errorStr);
            }
            
            if ( myFileManager.isFile() )
            {
                tmp.push_back( myFileManager.getFullFileName() );
            }
            else
            {
                myFileManager.setStringWithNamesOfFilesInDirectory( tmp );
            }
        }
        for (size_t i = 0; i < tmp.size(); i++)
        {
            vectorOfFileNames.push_back(tmp[i]);
        }
        std::set<std::string> s( vectorOfFileNames.begin(), vectorOfFileNames.end() );
        vectorOfFileNames.assign( s.begin(), s.end() );
    }
    else
    {
        // check that the file/path name has been correctly specified
        const std::string&  fn       = static_cast<const RlString&>( args[arg_index_files].getVariable()->getRevObject() ).getValue();
        RevBayesCore::RbFileManager myFileManager( fn );
        
        if ( !myFileManager.testFile() && !(myFileManager.testDirectory() && myFileManager.getFileName() == "" ) )
        {
            std::string errorStr = "Could not find filename: " + myFileManager.getFileName() + "\n";
            myFileManager.formatError(errorStr);
            throw RbException(errorStr);
        }
        
        // set up a vector of strings containing the name or names of the files to be read
        if ( myFileManager.isFile() )
        {
            vectorOfFileNames.push_back( myFileManager.getFullFileName() );
        }
        else
        {
            myFileManager.setStringWithNamesOfFilesInDirectory( vectorOfFileNames );
        }
    }
    
    TraceTree *rv;
    if ( treetype == "clock" )
    {
        if(nexus) rv = readTreesNexus(vectorOfFileNames, true, thin);
        else rv = readTrees(vectorOfFileNames, sep, true, thin);
    }
    else if ( treetype == "non-clock" )
    {
        RevBayesCore::Clade og;
        if ( args[arg_index_outgroup].getVariable() != NULL && args[arg_index_outgroup].getVariable()->getRevObject() != RevNullObject::getInstance())
        {
            og = static_cast<const Clade &>( args[arg_index_outgroup].getVariable()->getRevObject() ).getValue();
        }

        if(nexus) rv = readTreesNexus(vectorOfFileNames, false, thin);
        else rv = readTrees(vectorOfFileNames, sep, false, thin);
    }
    else
    {
        throw RbException("Unknown tree type to read.");
    }
    
    long burnin = 0;

    RevObject& b = args[arg_index_burnin].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burninFrac = static_cast<const Probability &>(b).getValue();
        burnin = long( floor( rv->getValue().size()*burninFrac ) );
    }

    rv->getValue().setBurnin(burnin);

    return new RevVariable( rv );
}


/** Get argument rules */
const ArgumentRules& Func_readTreeTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        //        std::vector<TypeSpec> branchRateTypeentRule( "branchRates", branchRateTypes, "The global or branch-specific rate multipliers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> fileTypes;
        fileTypes.push_back( RlString::getClassTypeSpec() );
        fileTypes.push_back( ModelVector<RlString>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "file"     , fileTypes, "The name of the tree trace file(s).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> tree_options;
        tree_options.push_back( "clock" );
        tree_options.push_back( "non-clock" );
        argumentRules.push_back( new OptionRule( "treetype", new RlString("clock"), tree_options, "The type of trees." ) );
        argumentRules.push_back( new ArgumentRule( "outgroup"   , Clade::getClassTypeSpec(), "The clade (consisting of one or more taxa) used as an outgroup.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "separator", RlString::getClassTypeSpec(), "The separator/delimiter between values in the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );

        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes     , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );

        argumentRules.push_back( new ArgumentRule( "thinning", Natural::getClassTypeSpec(), "The frequency of samples to read, i.e., we will only used every n-th sample where n is defined by this argument.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( 1l ) ) );

        argumentRules.push_back( new ArgumentRule( "nexus", RlBoolean::getClassTypeSpec(), "Whether the file to read is in NEXUS format.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readTreeTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_readTreeTrace";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readTreeTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readTreeTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTreeTrace";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readTreeTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readTreeTrace::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TraceTree::getClassTypeSpec();
    return return_typeSpec;
}


TraceTree* Func_readTreeTrace::readTrees(const std::vector<std::string> &vector_of_file_names, const std::string &delimitter, bool clock, long thinning)
{
    
    std::vector<RevBayesCore::TraceTree> data;
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.
    std::map<std::string,std::string> file_ap;
    for (std::vector<std::string>::const_iterator p = vector_of_file_names.begin(); p != vector_of_file_names.end(); ++p)
    {
        bool has_header_been_read = false;
        const std::string &fn = *p;
        
        RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(fn);
        
        // let us quickly count the number of lines
        size_t lines = 0;
        std::ifstream tmp_in_file( fm.getFullFileName().c_str() );
        
        if ( !tmp_in_file )
        {
            throw RbException( "Could not open file \"" + fn + "\"" );
        }
        while ( tmp_in_file.good() )
        {
            // Read a line
            std::string line;
            RevBayesCore::RbFileManager reader = RevBayesCore::RbFileManager();
            reader.safeGetline(tmp_in_file, line);
            if (line.length() == 0 || line[0] == '#')
            {
                continue;
            }
            ++lines;
        }
        tmp_in_file.close();
        
        RevBayesCore::ProgressBar progress = RevBayesCore::ProgressBar( lines, 0 );

        // now we actually process the input
        
        /* Open file */
        std::ifstream in_file( fm.getFullFileName().c_str() );
        
        if ( !in_file )
        {
            throw RbException( "Could not open file \"" + fn + "\"" );
        }
        
        /* Initialize */
        std::string commandLine;
        RBOUT( "Processing file \"" + fn + "\"");
        
        size_t n_samples = 0;
        size_t index = 0;
        progress.start();
        
        /* Command-processing loop */
        while ( in_file.good() )
        {
            
            // Read a line
            std::string line;
            RevBayesCore::RbFileManager reader = RevBayesCore::RbFileManager();
            reader.safeGetline(in_file, line);
            
            // skip empty lines
            //line = StringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
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
            
            // we should provide other delimiters too
            StringUtilities::stringSplit(line, delimitter, columns);
            
            // we assume a header at the first line of the file
            if ( has_header_been_read == false )
            {
                
                for (size_t j=1; j<columns.size(); j++)
                {
                    
                    std::string parmName = columns[j];
                    if ( parmName == "Posterior" || parmName == "Likelihood" || parmName == "Prior")
                    {
                        continue;
                    }
                    index = j;
                    
                    RevBayesCore::TraceTree t = RevBayesCore::TraceTree( clock );
                    
                    t.setParameterName(parmName);
                    t.setFileName(fn);
                    
                    data.push_back( t );
                }
                
                has_header_been_read = true;
                
                continue;
            }
            
            // increase our sample counter
            ++n_samples;
            
            // we need to check if we skip this sample in case of thinning.
            if ( (n_samples-1) % thinning > 0 )
            {
                continue;
            }
            
            // adding values to the Traces
            RevBayesCore::TraceTree& t = data[0];
            
            RevBayesCore::Tree *tau = NULL;
            if ( clock == true )
            {
                RevBayesCore::NewickConverter c;
                RevBayesCore::Tree *blTree = c.convertFromNewick( columns[index] );
                tau = RevBayesCore::TreeUtilities::convertTree( *blTree );
            }
            else
            {
                RevBayesCore::NewickConverter c;
                tau = c.convertFromNewick( columns[index] );
                
                // moved the reroot functionality below into
                // RevBayesCore::TreeSummary<BranchLengthTree>::summarizeTrees()
                // so that tree traces retain their original topology
                // will freyman 12/12/14
                
                //            if ( outgroup == "" )
                //            {
                //                RevBayesCore::BranchLengthTree& referenceTree = *tau;
                //                outgroup = referenceTree.getTipNode(0).getName();
                //            }
                //            // re-root the tree so that we can compare the the trees
                //            tau->reroot( outgroup );
                
                //                if ( outgroup.size() > 0 )
                //                {
                //                    //                tau->reroot( outgroup, false );
                //                }
            }
            
            t.addObject( tau );
            progress.update( n_samples );
            
        }
        in_file.close();
        
        progress.finish();
    }
    
    return new TraceTree( data[0] );
}

/** Create tree trace from one or several Nexus file(s)
 * @param fns vector of file names
 * @param clock whether trees have a clock
 * @param thin keep only every thin-th sample
 * @return tree trace
 * @see Func_readTrees::execute
 *
 * @note if multiple files are given, the traces will all be appended without regard for burnin
 * */
TraceTree *Func_readTreeTrace::readTreesNexus(const std::vector<string> &fns, bool clock, long thin) {
    RevBayesCore::TraceTree tt = RevBayesCore::TraceTree();
    tt.setParameterName("tree");

    for (size_t i=0; i<fns.size(); ++i)
    {
        const std::string fn = fns[i];
        // get the global instance of the NCL reader and clear warnings from its warnings buffer
        RevBayesCore::NclReader reader = RevBayesCore::NclReader();

        std::vector<RevBayesCore::Tree*> tmp;
        if ( clock ) {
            tmp = reader.readTimeTrees( fn );
        }
        else {
            tmp = *reader.readBranchLengthTrees( fn );
        }
        int nsamples = 0;
        for (size_t j=0; j<tmp.size(); ++j)
        {
            RevBayesCore::Tree* t = tmp[i];
            if(nsamples % thin == 0) tt.addObject(t);
            nsamples++;
        }
    }

    return new TraceTree(tt);
}


