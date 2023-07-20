#include "Func_collapseSA.h"

#include <stddef.h>

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
#include "GenericFunction.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TypeSpec.h"

using namespace RevLanguage;

RevBayesCore::Tree* collapseSampledAncestorsFunc(const RevBayesCore::Tree& tree)
{
    auto tree2 = tree.clone();

    tree2->collapseSampledAncestors();

    return tree2;
}

/** default constructor */
Func_collapseSA::Func_collapseSA( void ) : TypedFunction<Tree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_collapseSA* Func_collapseSA::clone( void ) const
{
    
    return new Func_collapseSA( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_collapseSA::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::Tree >* tree = static_cast<const Tree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::Tree >( collapseSampledAncestorsFunc, tree);
}


/* Get argument rules */
const ArgumentRules& Func_collapseSA::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_collapseSA::getClassType(void)
{
    
    static std::string rev_type = "Func_collapseSA";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_collapseSA::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_collapseSA::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCollapseSA";
    
    return f_name;
}


const TypeSpec& Func_collapseSA::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
