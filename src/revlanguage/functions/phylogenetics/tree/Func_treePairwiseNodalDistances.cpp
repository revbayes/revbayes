#include "Func_treePairwiseNodalDistances.h"

#include "DistanceMatrix.h"
#include "RlDeterministicNode.h"
#include "TreePairwiseNodalDistances.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlTree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;


/** default constructor */
Func_treePairwiseNodalDistances::Func_treePairwiseNodalDistances( void ) : TypedFunction< DistanceMatrix >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_treePairwiseNodalDistances * Func_treePairwiseNodalDistances::clone( void ) const
{
    
    return new Func_treePairwiseNodalDistances( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::DistanceMatrix >* Func_treePairwiseNodalDistances::createFunction( void ) const
{
    RevBayesCore::TypedDagNode < RevBayesCore::Tree >* tau = static_cast<const Tree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TreePairwiseNodalDistances * f = new RevBayesCore::TreePairwiseNodalDistances( tau );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_treePairwiseNodalDistances::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_treePairwiseNodalDistances::getClassType(void)
{
    
    static std::string rev_type = "Func_treePairwiseNodalDistances";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_treePairwiseNodalDistances::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<Tree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_treePairwiseNodalDistances::getFunctionName( void ) const
{
    std::string f_name = "fnTreePairwiseNodalDistances";
    
    return f_name;
}


const TypeSpec& Func_treePairwiseNodalDistances::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

