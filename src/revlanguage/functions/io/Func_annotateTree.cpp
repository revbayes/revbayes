#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_annotateTree.h"
#include "NexusWriter.h"
#include "Probability.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "RlTree.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_annotateTree* Func_annotateTree::clone( void ) const
{
    
    return new Func_annotateTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_annotateTree::execute( void )
{
    size_t arg_index = 0;
    
    // get the input tree
    RevBayesCore::Tree* tree = static_cast<const Tree&>( this->args[arg_index++].getVariable()->getRevObject() ).getValue().clone();
    
    // get the  tree trace(s)
    RevBayesCore::TreeSummary* summary = NULL;

    RevObject& b = args[arg_index++].getVariable()->getRevObject();
    if ( b.isType( TraceTree::getClassTypeSpec() ) )
    {
        summary = &( static_cast<TraceTree &>(b).getValue() );
    }
    else
    {
        std::vector<TraceTree>& tt = static_cast<const WorkspaceVector<TraceTree> &>(b).getValue();

        std::vector<RevBayesCore::TraceTree* > vec;
        for(size_t i = 0; i < tt.size(); i++)
        {
            vec.push_back(&(tt[i].getValue()));
        }

        summary = new RevBayesCore::TreeSummary( vec, tt.front().getValue().isClock() );
    }
    
    // get the filename
    const std::string& filename = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    

    RevBayesCore::TraceTree::AnnotationReport report;

    report.node_ages               = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.conditional_clade_ages  = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.conditional_clade_probs = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.conditional_tree_ages   = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.node_ages_HPD           = static_cast<const Probability &>(args[arg_index++].getVariable()->getRevObject()).getValue();
    report.mean_node_ages          = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    report.sampled_ancestor_probs  = static_cast<const RlBoolean &>( this->args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    // do not make a new tree summary object
    // this way we don't need to resummarize every time we annotate a tree
    // RevBayesCore::TreeSummary summary = RevBayesCore::TreeSummary( tt.getValue() );
    
    bool verbose = true;
    summary->annotateTree( *tree, report, verbose );
    
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
    
    Tree* t;
    if ( summary->isClock() )
    {
        t = new TimeTree( tree );
    }
    else
    {
        t = new BranchLengthTree( tree );
    }

    if ( b.isType( TraceTree::getClassTypeSpec() ) == false )
    {
        delete summary;
    }

    return new RevVariable( t );
}



/** Get argument rules */
const ArgumentRules& Func_annotateTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec()        , "The input tree which will be annotated.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        std::vector<TypeSpec> traceArgRules;
        traceArgRules.push_back( TraceTree::getClassTypeSpec() );
        traceArgRules.push_back( WorkspaceVector<TraceTree>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "trace", traceArgRules, "The sample trace.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec()    , "The name of the file where to store the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("")  ) );

        argumentRules.push_back( new ArgumentRule( "ages" , RlBoolean::getClassTypeSpec() , "Annotate node ages?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "ccAges" , RlBoolean::getClassTypeSpec() , "Annotate conditional clade ages?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "ccp" , RlBoolean::getClassTypeSpec() , "Annotate conditional clade probabilities?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "conditionalAges" , RlBoolean::getClassTypeSpec() , "Annotate node ages conditional on the input topology?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "hpd"   ,    Probability::getClassTypeSpec() , "The probability mass of the highest posterior density node age interval.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95) ) );
        argumentRules.push_back( new ArgumentRule( "mean" , RlBoolean::getClassTypeSpec() , "Annotate node ages using the mean age instead of the median?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "sampledAncestors" , RlBoolean::getClassTypeSpec() , "Annotate sampled ancestor probs?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_annotateTree::getClassType(void)
{
    
    static std::string rev_type = "Func_annotateTree<" + Tree::getClassType() + ">";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_annotateTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_annotateTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "annotateTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_annotateTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_annotateTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Tree::getClassTypeSpec();
    return return_typeSpec;
}

