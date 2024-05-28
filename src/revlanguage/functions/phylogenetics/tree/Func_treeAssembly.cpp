#include "Func_treeAssembly.h"

#include <cstddef>

#include "ModelVector.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlTree.h"
#include "RlDeterministicNode.h"
#include "TreeAssemblyFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_treeAssembly::Func_treeAssembly( void ) : TypedFunction<Tree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_treeAssembly* Func_treeAssembly::clone( void ) const
{
    
    return new Func_treeAssembly( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_treeAssembly::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // TreeAssemblyFunction acts directly on the value of the topology node
    // Current topology variable cannot already be parent of another TreeAssemblyFunction
    bool topology_in_use = false;
    const std::vector<RevBayesCore::DagNode*>& tauChildren = tau->getChildren();
    for (size_t i = 0; i < tauChildren.size(); i++)
    {
        RevBayesCore::DeterministicNode<RevBayesCore::Tree>* tauChild = dynamic_cast<RevBayesCore::DeterministicNode<RevBayesCore::Tree>*>(tauChildren[i]);
        if ( tauChild != NULL )
        {
            RevBayesCore::TreeAssemblyFunction* tf = dynamic_cast<RevBayesCore::TreeAssemblyFunction*>(&tauChild->getFunction());
            if ( tf != NULL )
            {
                topology_in_use = true;
            }
        }
    }
    if ( topology_in_use == true )
    {
        throw RbException("Variable \"" + tau->getName() + "\" cannot be used with more than one treeAssembly function.");
    }
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* brlens = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::RbVector<double> tmp = brlens->getValue();
    RevBayesCore::TreeAssemblyFunction* f = new RevBayesCore::TreeAssemblyFunction( tau, brlens );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_treeAssembly::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "topology", Tree::getClassTypeSpec(), "The tree topology variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "brlens",   ModelVector<RealPos>::getClassTypeSpec(), "The vector of branch lengths.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_treeAssembly::getClassType(void)
{
    
    static std::string rev_type = "Func_treeAssembly";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_treeAssembly::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_treeAssembly::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "treeAssembly" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_treeAssembly::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTreeAssembly";
    
    return f_name;
}


const TypeSpec& Func_treeAssembly::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
