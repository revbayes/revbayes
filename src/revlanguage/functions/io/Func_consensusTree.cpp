#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Func_consensusTree.h"
#include "ModelVector.h"
#include "NexusWriter.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "RlUtils.h"
#include "Probability.h"
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
Func_consensusTree* Func_consensusTree::clone(void) const
{
    
    return new Func_consensusTree(*this);
}


/** Execute function */
RevPtr<RevVariable> Func_consensusTree::execute(void)
{
    size_t arg_index = 0;

    TraceTree& tt = static_cast<TraceTree&>( args[arg_index++].getVariable()->getRevObject() );
    
    double cutoff = static_cast<const Probability &>(args[arg_index++].getVariable()->getRevObject()).getValue();

    const std::string& filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();


    RevBayesCore::TraceTree::AnnotationReport report;

    report.node_ages_HPD                 = static_cast<const Probability &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    report.mean_node_ages                = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.sampled_ancestor_probs        = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.force_positive_branch_lengths = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();

    bool verbose = true;
    RevBayesCore::Tree* tree = tt.getValue().mrTree(report, cutoff, verbose);
    
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
const ArgumentRules& Func_consensusTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "trace", TraceTree::getClassTypeSpec(), "The trace of tree samples.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "cutoff"   , Probability::getClassTypeSpec()  , "The minimum threshold for clade probabilities.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.5) ) );
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec() , "The name of the file for storing the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        
        argumentRules.push_back( new ArgumentRule( "hpd"   ,    Probability::getClassTypeSpec() , "The probability mass of the highest posterior density node age interval.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95) ) );
        argumentRules.push_back( new ArgumentRule( "mean" , RlBoolean::getClassTypeSpec() , "Annotate node ages using the mean age instead of the median?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "sampledAncestors" , RlBoolean::getClassTypeSpec() , "Annotate sampled ancestor probs?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "positiveBranchLengths" , RlBoolean::getClassTypeSpec() , "Force negative branch lengths to be short but positive?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_consensusTree::getClassType(void)
{
    
    static std::string rev_type = "Func_consensusTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_consensusTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_consensusTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "consensusTree";
    
    return f_name;
}

/**
 * Get the name for this procedure.
 */
std::vector<std::string> Func_consensusTree::getFunctionNameAliases( void ) const
{
    std::vector<std::string> aliases;

    aliases.push_back("conTree");
    aliases.push_back("sumt");

    return aliases;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_consensusTree::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Seraina Klopfstein" );
    authors.push_back( "Will Freyman" );
    authors.push_back( "Will Pett" );
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_consensusTree::getHelpDescription(void) const
{
    std::string description = "";
    description += "Calculates the majority-rule consensus topology from a trace of trees and summarizes branch lengths.";
    
    return description;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_consensusTree::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Read in tree trace\n";
    example += "tree_trace = readTreeTrace(\"output/my.trees\", burnin=0.25)\n";
    example += "\n";
    example += "# Generate the majority-rule consensus tree\n";
    example += "map_tree = consensusTree(trace=tree_trace, cutoff=0.5, file=\"consensus.tree\")\n";
    
    return example;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_consensusTree::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mapTree" );
    see_also.push_back( "mccTree" );
    see_also.push_back( "treeTrace" );
    see_also.push_back( "readTreeTrace" );
    
    
    return see_also;
}


/** Get type spec */
const TypeSpec& Func_consensusTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_consensusTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    
    return return_typeSpec;
}

