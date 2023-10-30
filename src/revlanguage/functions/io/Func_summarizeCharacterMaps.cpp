#include "Func_summarizeCharacterMaps.h"

#include <cmath>
#include <cstddef>
#include <sstream>
#include <vector>

#include "ArgumentRule.h"
#include "Delimiter.h"
#include "JointAncestralStateTrace.h"
#include "Probability.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlTraceTree.h"
#include "RlAncestralStateTrace.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Integer.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "RlTree.h"
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
Func_summarizeCharacterMaps* Func_summarizeCharacterMaps::clone( void ) const
{
    
    return new Func_summarizeCharacterMaps( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_summarizeCharacterMaps::execute( void )
{
    
    // get the input summary tree
    RevBayesCore::Tree input_tree = RevBayesCore::Tree();
    if (this->args[0].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        const RevBayesCore::TypedDagNode<RevBayesCore::Tree> *input_tree_dag = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
        input_tree = input_tree_dag->getValue();
    }

    // get the vector of stochastic character map traces
    const WorkspaceVector<AncestralStateTrace>& ast_vector = static_cast<const WorkspaceVector<AncestralStateTrace> &>( args[1].getVariable()->getRevObject() );
    std::vector<RevBayesCore::AncestralStateTrace> ancestralstate_traces;
    for (int i = 0; i < ast_vector.size(); ++i)
    {
        ancestralstate_traces.push_back( ast_vector[i].getValue() );
    }
    
    // get the tree trace
    const TraceTree& tt = static_cast<const TraceTree&>( args[2].getVariable()->getRevObject() );
    
    // make a new tree summary object, using the tree trace if necessary
    RevBayesCore::TraceTree tree_trace;
    if (args[2].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        tree_trace = tt.getValue();
    }
    
    // get the filename to write output
    const std::string& filename = static_cast<const RlString&>( args[3].getVariable()->getRevObject() ).getValue();
   
    // check if burnin was entered as integer or probability
    int burnin = 0;
    RevObject& b = args[4].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = (int)static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burninFrac = static_cast<const Probability &>(b).getValue();
        burnin = int( floor( ancestralstate_traces[0].size() * burninFrac ) );
    }
    
    const std::string& sep = static_cast<const RlString  &>( args[5].getVariable()->getRevObject() ).getValue();

    bool verbose = static_cast<const RlBoolean &>( args[6].getVariable()->getRevObject() ).getValue();
    
    bool bwd_time = static_cast<const RlBoolean &>( args[7].getVariable()->getRevObject() ).getValue();
    
    // summarize stochastic character maps
    RevBayesCore::JointAncestralStateTrace joint_trace(ancestralstate_traces, tree_trace);
    joint_trace.setBurnin(burnin);

    joint_trace.summarizeCharacterMaps(input_tree, filename, verbose, sep, bwd_time);

    return NULL;
}



/** Get argument rules */
const ArgumentRules& Func_summarizeCharacterMaps::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The input tree to summarize ancestral states over.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "character_map_trace_vector", WorkspaceVector<AncestralStateTrace>::getClassTypeSpec(), "A vector of stochastic character map traces.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "tree_trace", TraceTree::getClassTypeSpec(), "A trace of tree samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec() , "The name of the file to store the summarized character histories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes  , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        argumentRules.push_back( new Delimiter() );
        argumentRules.push_back( new ArgumentRule( "verbose"   , RlBoolean::getClassTypeSpec()  , "Printing verbose output", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "backward_time"   , RlBoolean::getClassTypeSpec()  , "Printing verbose output", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_summarizeCharacterMaps::getClassType(void)
{
    
    static std::string rev_type = "Func_summarizeCharacterMaps";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_summarizeCharacterMaps::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_summarizeCharacterMaps::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "summarizeCharacterMaps";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_summarizeCharacterMaps::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_summarizeCharacterMaps::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}

