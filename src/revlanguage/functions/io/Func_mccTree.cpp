#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Func_mccTree.h"
#include "ModelVector.h"
#include "NexusWriter.h"
#include "Probability.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "RlUtils.h"
#include "StringUtilities.h"

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
Func_mccTree* Func_mccTree::clone( void ) const
{
    
    return new Func_mccTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_mccTree::execute( void )
{
    size_t arg_index = 0;

    TraceTree& tt = static_cast<TraceTree&>( args[arg_index++].getVariable()->getRevObject() );

    const std::string& filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    
    RevBayesCore::TraceTree::AnnotationReport report;

    report.conditional_clade_ages        = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.conditional_clade_probs       = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.conditional_tree_ages         = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.node_ages_HPD                 = static_cast<const Probability &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    report.mean_node_ages                = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.sampled_ancestor_probs        = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.force_positive_branch_lengths = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    bool verbose = true;
    RevBayesCore::Tree* tree = tt.getValue().mccTree(report, verbose);
    
    
    if ( filename != "" )
    {
        
        RevBayesCore::NexusWriter writer(filename);
        writer.openStream(false);
        
        std::vector<RevBayesCore::Taxon> taxa;
        tree->getRoot().getTaxa(taxa);
        RevBayesCore::Clade c( taxa );
        writer.writeNexusBlock(c);
        
        writer.writeNexusBlock(*tree);
        
        writer.closeStream();
        
    }
    
    Tree* t;
    if ( tt.getValue().isClock() )
    {
        t = new TimeTree( tree );
    }
    else
    {
        t = new BranchLengthTree( tree );
    }

    return new RevVariable( t );
}



/** Get argument rules */
const ArgumentRules& Func_mccTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "trace", TraceTree::getClassTypeSpec(), "The samples of trees from the posterior.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "The name of the file where to store the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("")  ) );
        
        argumentRules.push_back( new ArgumentRule( "ccAges" , RlBoolean::getClassTypeSpec() , "Annotate conditional clade ages?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "ccp" , RlBoolean::getClassTypeSpec() , "Annotate conditional clade probabilities?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "conditionalAges" , RlBoolean::getClassTypeSpec() , "Annotate node ages conditional on the topology?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "hpd"   ,    Probability::getClassTypeSpec() , "The probability mass of the highest posterior density node age interval.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95) ) );
        argumentRules.push_back( new ArgumentRule( "mean" , RlBoolean::getClassTypeSpec() , "Annotate node ages using the mean age instead of the median?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "sampledAncestors" , RlBoolean::getClassTypeSpec() , "Annotate sampled ancestor probs?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "positiveBranchLengths" , RlBoolean::getClassTypeSpec() , "Force negative branch lengths to be short but positive?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_mccTree::getClassType(void)
{
    
    static std::string revType = "Func_mccTree";
    
    return revType;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_mccTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_mccTree::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Will Pett" );
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_mccTree::getHelpDescription(void) const
{
    std::string description = "";
    description += "Finds the maximum clade credibility (MCC) topology from a trace of trees and summarizes branch lengths.";
    
    return description;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_mccTree::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Read in tree trace\n";
    example += "tree_trace = readTreeTrace(\"output/my.trees\", burnin=0.25)\n";
    example += "\n";
    example += "# Generate the MCC tree\n";
    example += "map_tree = mccTree(trace=tree_trace, file=\"mcc.tree\")\n";
    
    return example;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_mccTree::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "consensusTree" );
    see_also.push_back( "mapTree" );
    see_also.push_back( "treeTrace" );
    see_also.push_back( "readTreeTrace" );
    
    
    return see_also;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_mccTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mccTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_mccTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_mccTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}

