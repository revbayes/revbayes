#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "OptionRule.h"
#include "Func_ancestralStateTree.h"
#include "JointAncestralStateTrace.h"
#include "NexusWriter.h"
#include "Probability.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlTraceTree.h"
#include "RlAncestralStateTrace.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Integer.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "RlTree.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Trace.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_ancestralStateTree* Func_ancestralStateTree::clone( void ) const
{
    
    return new Func_ancestralStateTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_ancestralStateTree::execute( void )
{
    
    // get the input tree
    const RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tree_node = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // get vector of ancestral state traces
    const WorkspaceVector<AncestralStateTrace>& ast_vector = static_cast<const WorkspaceVector<AncestralStateTrace> &>( args[1].getVariable()->getRevObject() );
    std::vector<RevBayesCore::AncestralStateTrace> ancestralstate_traces;
    for (int i = 0; i < ast_vector.size(); ++i)
    {
        ancestralstate_traces.push_back( ast_vector[i].getValue() );
    }
    
    // get the ancestral state tree trace
    const TraceTree& tt = static_cast<const TraceTree&>( args[2].getVariable()->getRevObject() );
    
    // make a new tree summary object
    RevBayesCore::TraceTree tree_trace;

    if (args[2].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        tree_trace = tt.getValue();
    }
    
    // should we annotate start states?
    bool start_states = static_cast<const RlBoolean &>(args[3].getVariable()->getRevObject()).getValue();
    
    // get the filename
    const std::string& filename = static_cast<const RlString&>( args[4].getVariable()->getRevObject() ).getValue();
    
    int burnin = 0;

    RevObject& b = args[5].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = (int)static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burninFrac = static_cast<const Probability &>(b).getValue();
        burnin = int( floor( ancestralstate_traces[0].size() * burninFrac ) );
    }

    std::string summary_stat = static_cast<const RlString&>( args[6].getVariable()->getRevObject() ).getValue();
    

    std::string reconstruction = static_cast<const RlString &>(args[7].getVariable()->getRevObject()).getValue();
    bool conditional = false;
    if ( reconstruction == "conditional" )
    {
        conditional = true;
    }
    if ( reconstruction == "joint" )
    {
        throw RbException("Joint ancestral state summaries are not yet implemented. Coming soon!");
    }
    
    int site = (int)static_cast<const Integer &>(args[8].getVariable()->getRevObject()).getValue();
    if ( site == 0 )
    {
        throw RbException("In Rev we index using a base '1'. That means, the first site has position '1'. You entered '0' for the site index.");
    }
    --site;
    
    bool verbose = static_cast<const RlBoolean &>(args[9].getVariable()->getRevObject()).getValue();
    
    // get the tree with ancestral states
    RevBayesCore::JointAncestralStateTrace joint_trace(ancestralstate_traces, tree_trace);
    joint_trace.setBurnin(burnin);

    RevBayesCore::Tree* tree;
    if (start_states)
    {
        tree = joint_trace.cladoAncestralStateTree(tree_node->getValue(), summary_stat, site, conditional, false, verbose);
    }
    else
    {
        size_t num_states = 3;
        num_states = static_cast<const Natural &>(args[10].getVariable()->getRevObject()).getValue();

        tree = joint_trace.ancestralStateTree(tree_node->getValue(), summary_stat, num_states, site, conditional, false, verbose);
    }
    
    // return the tree
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
    
    return new RevVariable( new Tree( tree ) );
}



/** Get argument rules */
const ArgumentRules& Func_ancestralStateTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The input tree to summarize ancestral states over.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "ancestral_state_trace_vector", WorkspaceVector<AncestralStateTrace>::getClassTypeSpec(), "A vector of ancestral state traces.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "tree_trace", TraceTree::getClassTypeSpec(), "A trace of tree samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "include_start_states", RlBoolean::getClassTypeSpec(), "Annotate start states as well as end states for each branch. Only applicable for cladogenetic processes.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec() , "The name of the file to store the annotated tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes  , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        std::vector<std::string> summary_stats;
        summary_stats.push_back( "MAP" );
        summary_stats.push_back( "mean" );
        argumentRules.push_back( new OptionRule( "summary_statistic", new RlString("MAP"), summary_stats, "The statistic used to summarize ancestral states. 'MAP' displays the 3 states with highest posterior probabilities. 'mean' displays the mean value and 95% CI." ) );

        std::vector<std::string> reconstruction;
        reconstruction.push_back( "conditional" );
        reconstruction.push_back( "joint" );
        reconstruction.push_back( "marginal" );
        argumentRules.push_back( new OptionRule( "reconstruction", new RlString("marginal"), reconstruction, "'joint' and 'conditional' should only be used to summarize ancestral states sampled from the joint distribution. 'marginal' can be used for states sampled from the joint or marginal distribution." ) );
        argumentRules.push_back( new ArgumentRule( "site"     , Integer::getClassTypeSpec()  , "The character site to be summarized.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Integer(1) ) );
        argumentRules.push_back( new ArgumentRule( "verbose"   , RlBoolean::getClassTypeSpec()  , "Printing verbose output", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "nStates"   , Natural::getClassTypeSpec()  , "The number of states for which we compute the posterior probability. By default it will be the best three.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(3) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_ancestralStateTree::getClassType(void)
{
    
    static std::string rev_type = "Func_ancestralStateTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_ancestralStateTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_ancestralStateTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "ancestralStateTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_ancestralStateTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_ancestralStateTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}

