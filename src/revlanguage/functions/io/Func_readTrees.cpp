#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Func_readTrees.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "NewickConverter.h"
#include "OptionRule.h"
#include "RbFileManager.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "Tree.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readTrees* Func_readTrees::clone( void ) const
{

    return new Func_readTrees( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readTrees::execute( void )
{

    // get the information from the arguments for reading the file
    const std::string&  fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  text = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string&  treetype = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    bool unroot_nonclock = static_cast<const RlBoolean&>( args[3].getVariable()->getRevObject() ).getValue();

    if (fn != "")
    {
        // get the global instance of the NCL reader and clear warnings from its warnings buffer
        RevBayesCore::NclReader reader = RevBayesCore::NclReader();

        if ( treetype == "clock" )
        {
            ModelVector<TimeTree> *trees = new ModelVector<TimeTree>();
            
            std::vector<RevBayesCore::Tree*>* tmp = reader.readTimeTrees( fn );
            for (std::vector<RevBayesCore::Tree*>::iterator t = tmp->begin(); t != tmp->end(); ++t)
            {
                trees->push_back( TimeTree(*t) );
            }
            delete tmp;
            
            return new RevVariable( trees );
        }
        else if ( treetype == "non-clock" )
        {
            ModelVector<BranchLengthTree> *trees = new ModelVector<BranchLengthTree>();
            
            std::vector<RevBayesCore::Tree*>* tmp_trees = reader.readBranchLengthTrees( fn );
            // FIXME - this can be NULL for branch-length trees, but not for time trees!  Bad.
            if (tmp_trees)
            {
                for (auto tree: *tmp_trees)
                {
                    if (unroot_nonclock)
                    {
                        tree->removeRootIfDegree2();
//                      Perhaps we should mark the tree unrooted, since we have removed the old root,
//                        and chosen a neighbor as the now root.
//                      However, RevBayes has bugs with unrooted trees and may crash.
//                        tree->setRooted(false);
                    }

                    trees->push_back( BranchLengthTree(*tree) );
                }
                delete tmp_trees;
            }
            return new RevVariable( trees );
        }
        
    }

    if (text != "")
    {

        string aux;
        std::istringstream iss(text);
        RevBayesCore::NewickConverter c;
        
        if ( treetype == "clock" )
        {
            ModelVector<TimeTree> *trees = new ModelVector<TimeTree>();
            while (RevBayesCore::safeGetline(iss, aux))
            {
                RevBayesCore::Tree *blTree = c.convertFromNewick( aux );
                trees->push_back( TimeTree(*blTree) );
                
                delete blTree;
            }
            return new RevVariable( trees );
        }
        else if ( treetype == "non-clock" )
        {
            ModelVector<BranchLengthTree> *trees = new ModelVector<BranchLengthTree>();
            while (RevBayesCore::safeGetline(iss, aux))
            {
                RevBayesCore::Tree *blTree = c.convertFromNewick( aux );

                if (unroot_nonclock)
                {
                    blTree->removeRootIfDegree2();
                    blTree->setRooted(false);
                }

                trees->push_back( BranchLengthTree(*blTree) );
                
                delete blTree;
            }
            return new RevVariable( trees );
            
        }

        
        
    }


    return NULL;

}


/** Get argument rules */
const ArgumentRules& Func_readTrees::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {

        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file containing the trees.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        argumentRules.push_back( new ArgumentRule( "text", RlString::getClassTypeSpec(), "A string containing one or several newick trees, separated by end of lines.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        
        std::vector<std::string> tree_options;
        tree_options.push_back( "clock" );
        tree_options.push_back( "non-clock" );
        argumentRules.push_back( new OptionRule( "treetype", new RlString("clock"), tree_options, "The type of trees." ) );
        
        argumentRules.push_back( new ArgumentRule( "unroot_nonclock", RlBoolean::getClassTypeSpec(), "Remove a degree-2 root node and set the tree unrooted, if treetype is non-clock.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        rules_set = true;
    }

    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readTrees::getClassType(void)
{

    static std::string rev_type = "Func_readTrees";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readTrees::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readTrees::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTrees";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readTrees::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readTrees::getReturnType( void ) const
{

    static TypeSpec return_typeSpec = ModelVector<TimeTree>::getClassTypeSpec();
    return return_typeSpec;
}
