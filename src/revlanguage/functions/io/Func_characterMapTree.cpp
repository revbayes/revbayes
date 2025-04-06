#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "OptionRule.h"
#include "Func_characterMapTree.h"
#include "JointAncestralStateTrace.h"
#include "Probability.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlTraceTree.h"
#include "RlAncestralStateTrace.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Integer.h"
#include "RbBoolean.h"
#include "RbFileManager.h"
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
Func_characterMapTree* Func_characterMapTree::clone( void ) const
{
    
    return new Func_characterMapTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_characterMapTree::execute( void )
{
    
    size_t arg_index = 0;
    
    // get the input tree
    const RevBayesCore::TypedDagNode<RevBayesCore::Tree> *it = static_cast<const Tree&>( this->args[arg_index++].getVariable()->getRevObject() ).getDagNode();
    
    // get vector of ancestral state traces
    const WorkspaceVector<AncestralStateTrace>& ast_vector = static_cast<const WorkspaceVector<AncestralStateTrace> &>( args[arg_index++].getVariable()->getRevObject() );
    std::vector<RevBayesCore::AncestralStateTrace> ancestralstate_traces;
    for (int i = 0; i < ast_vector.size(); ++i)
    {
        ancestralstate_traces.push_back( ast_vector[i].getValue() );
    }
    
    // make a new tree summary object
    RevBayesCore::TraceTree tree_trace;
    const Argument& tt_arg = args[arg_index++];
    if (tt_arg.getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        
        // get the ancestral state tree trace
        const TraceTree& tt = static_cast<const TraceTree&>( tt_arg.getVariable()->getRevObject() );
        tree_trace = tt.getValue();
    }
    
    // get the filename for the tree with MAP character history
    RevBayesCore::path map_filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    // get the filename for the tree with posteriors for the MAP character history
    RevBayesCore::path map_pp_filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();

    // get the filename for the tree with shift probability for character history
    RevBayesCore::path map_shift_pp_filename = "";
    const Argument& map_shift_pp_filename_arg = args[arg_index++];
    if ( map_shift_pp_filename_arg.getVariable()->getRevObject() != RevNullObject::getInstance() )
    {
        map_shift_pp_filename = static_cast<const RlString&>( map_shift_pp_filename_arg.getVariable()->getRevObject() ).getValue();
    }
    
    int burnin;
    RevObject& b = args[arg_index++].getVariable()->getRevObject();
    if ( b.isType( Integer::getClassTypeSpec() ) )
    {
        burnin = (int)static_cast<const Integer &>(b).getValue();
    }
    else
    {
        double burnin_frac = static_cast<const Probability &>(b).getValue();
        burnin = int( floor( ancestralstate_traces[0].size() * burnin_frac ) );
    }

    std::string reconstruction = static_cast<const RlString &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    bool conditional = false;
    if ( reconstruction == "conditional" )
    {
        conditional = true;
    }
    if ( reconstruction == "joint" )
    {
        throw RbException("Joint ancestral state summaries are not yet implemented. Coming soon!");
    }
    
    int num_time_slices = (int)static_cast<const Integer &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    
    bool verbose = static_cast<const RlBoolean &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    
    // get the tree with ancestral states
    RevBayesCore::JointAncestralStateTrace joint_trace(ancestralstate_traces, tree_trace);
    joint_trace.setBurnin(burnin);
    RevBayesCore::Tree* tree = joint_trace.characterMapTree(it->getValue(), num_time_slices, conditional, false, verbose);
    
    // write the SIMMAP newick strings
    std::ofstream out_stream;
    if ( not map_filename.empty() )
    {
        RevBayesCore::createDirectoryForFile( map_filename );
        
        RevBayesCore::Tree t = RevBayesCore::Tree(*tree);
        t.renameNodeParameter("map_character_history", "character_history");
        out_stream.open( map_filename.string(), std::fstream::out);
        out_stream << t.getSimmapNewickRepresentation();
        out_stream.close();
    }
    if ( not map_pp_filename.empty() )
    {
        RevBayesCore::createDirectoryForFile( map_pp_filename );
        
        RevBayesCore::Tree t = RevBayesCore::Tree(*tree);
        t.renameNodeParameter("map_character_history_posteriors", "character_history");
        out_stream.open( map_pp_filename.string(), std::fstream::out);
        out_stream << t.getSimmapNewickRepresentation();
        out_stream.close();
    }
    if ( not map_shift_pp_filename.empty() )
    {
        RevBayesCore::createDirectoryForFile( map_shift_pp_filename );
        
        RevBayesCore::Tree t = RevBayesCore::Tree(*tree);
        t.renameNodeParameter("map_character_history_shift_prob", "character_history");
        out_stream.open( map_shift_pp_filename.string(), std::fstream::out);
        out_stream << t.getSimmapNewickRepresentation();
        out_stream.close();
    }
    
    // return the tree with annotations
    return new RevVariable( new Tree( tree ) );
}



/** Get argument rules */
const ArgumentRules& Func_characterMapTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        

        argumentRules.push_back( new ArgumentRule( "tree",                         Tree::getClassTypeSpec(),                                 "The input tree to summarize the character history over.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "ancestral_state_trace_vector", WorkspaceVector<AncestralStateTrace>::getClassTypeSpec(), "A vector of ancestral state traces.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "tree_trace",                   TraceTree::getClassTypeSpec(),                            "A trace of tree samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "character_file",               RlString::getClassTypeSpec(),                             "The name of the file to store the tree annotated with the MAP character history.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "posterior_file",               RlString::getClassTypeSpec(),                             "The name of the file to store the tree annotated with the posterior probabilities for the MAP character history.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "shift_prob_file",              RlString::getClassTypeSpec(),                             "The name of the file to store the tree annotated with the shift probabilities of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

		std::vector<TypeSpec> burnin_types;
        burnin_types.push_back( Probability::getClassTypeSpec() );
        burnin_types.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burnin_types  , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
   
        std::vector<std::string> reconstruction;
        reconstruction.push_back( "conditional" );
        reconstruction.push_back( "joint" );
        reconstruction.push_back( "marginal" );
        argumentRules.push_back( new OptionRule( "reconstruction", new RlString("marginal"), reconstruction, "'joint' and 'conditional' should only be used to summarize character maps sampled from the joint distribution. 'marginal' can be used for character maps sampled from the joint or marginal distribution." ) );

        argumentRules.push_back( new ArgumentRule( "num_time_slices",              Integer::getClassTypeSpec(),                              "The number of time slices to discretize the character history. Should be the same as used for the numeric ODE.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Integer(500) ) );
        argumentRules.push_back( new ArgumentRule( "verbose",                      RlBoolean::getClassTypeSpec(),                            "Printing verbose output", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_characterMapTree::getClassType(void)
{
    
    static std::string rev_type = "Func_characterMapTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_characterMapTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_characterMapTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "characterMapTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_characterMapTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_characterMapTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}

