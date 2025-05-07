#include "Func_pruneTree.h"

#include <cstddef>

#include "ModelVector.h"
#include "PruneTreeFunction.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "RlDeterministicNode.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "RlTree.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_pruneTree::Func_pruneTree( void ) : TypedFunction<Tree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_pruneTree* Func_pruneTree::clone( void ) const
{
    
    return new Func_pruneTree( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_pruneTree::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree&>( args[0].getVariable()->getRevObject() ).getDagNode();

    std::vector<RevBayesCore::Taxon> taxa;

    if( args[1].getVariable()->getRevObject().isType( ModelVector<Taxon>::getClassTypeSpec() ) )
    {
        taxa = static_cast<const ModelVector<Taxon>&>( args[1].getVariable()->getRevObject() ).getValue();
    }
    else
    {
        std::vector<std::string> strings = static_cast<const ModelVector<RlString>&>( args[1].getVariable()->getRevObject() ).getValue();

        for(size_t i = 0; i < strings.size(); i++)
        {
            taxa.push_back(RevBayesCore::Taxon(strings[i]));
        }
    }

    bool prune_fossils = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();

    RevBayesCore::PruneTreeFunction* f = new RevBayesCore::PruneTreeFunction( tau, taxa, args[1].getLabel() == "retain", prune_fossils );
    
    return f;
}


/** Execute function */
RevPtr<RevVariable> Func_pruneTree::execute( void )
{
    RevBayesCore::TypedFunction<RevBayesCore::Tree>* d = createFunction();
    RevBayesCore::DeterministicNode<RevBayesCore::Tree>* rv;
    rv = new DeterministicNode<RevBayesCore::Tree>("", d, this->clone());

    RevObject& ro = args[0].getVariable()->getRevObject();
    if ( ro.isType( TimeTree::getClassTypeSpec() ) )
    {
        return new RevVariable( new TimeTree(rv) );
    }
    else if ( ro.isType( BranchLengthTree::getClassTypeSpec() ) )
    {
        return new RevVariable( new BranchLengthTree(rv) );
    }

    return new RevVariable( new Tree(rv) );
}


/* Get argument rules */
const ArgumentRules& Func_pruneTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<std::string> labels;
        labels.push_back("prune");
        labels.push_back("retain");

        std::vector<TypeSpec> types;
        types.push_back(ModelVector<Taxon>::getClassTypeSpec());
        types.push_back(ModelVector<RlString>::getClassTypeSpec());

        argumentRules.push_back( new ArgumentRule( labels , types , "Taxon set to prune/retain in the tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "pruneFossils" , RlBoolean::getClassTypeSpec() , "Prune all fossils from tree?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_pruneTree::getClassType(void)
{
    
    static std::string rev_type = "Func_pruneTree";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_pruneTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_pruneTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPruneTree";
    
    return f_name;
}


const TypeSpec& Func_pruneTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
