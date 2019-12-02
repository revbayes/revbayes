#include "Func_treePairwiseDistances.h"
#include "DistanceMatrix.h"
#include "ModelVector.h"
#include "RlTimeTree.h"
#include "RlBranchLengthTree.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TreePairwiseDistances.h"
#include "TypedDagNode.h"

using namespace RevLanguage;


/** default constructor */
Func_treePairwiseDistances::Func_treePairwiseDistances( void ) : TypedFunction< DistanceMatrix >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_treePairwiseDistances * Func_treePairwiseDistances::clone( void ) const
{
    
    return new Func_treePairwiseDistances( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::DistanceMatrix >* Func_treePairwiseDistances::createFunction( void ) const
{
    
    // RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode < RevBayesCore::Tree >* tau = static_cast<const Tree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TreePairwiseDistances * f = new RevBayesCore::TreePairwiseDistances( tau );
    
    /*	if ( this->args[0].getVariable()->getRequiredTypeSpec().isDerivedOf( TimeTree::getClassTypeSpec() ) )
     {
     RevBayesCore::TypedDagNode < RevBayesCore::Tree >* tau = static_cast<const treeType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
     f = new RevBayesCore::TreePairwiseDistances < RevBayesCore::Tree > ( tau );
     }
     else if ( this->args[0].getVariable()->getRequiredTypeSpec().isDerivedOf( BranchLengthTree::getClassTypeSpec() ) )
     {
     RevBayesCore::TypedDagNode< RevBayesCore::BranchLengthTree >* tau = static_cast<const BranchLengthTree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
     f = new RevBayesCore::TreePairwiseDistances< RevBayesCore::BranchLengthTree > ( tau );
     }*/
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_treePairwiseDistances::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        /*std::vector<TypeSpec> treeTypes;
         treeTypes.push_back( TimeTree::getClassTypeSpec() );
         treeTypes.push_back( BranchLengthTree::getClassTypeSpec() );*/
        
        
        argumentRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_treePairwiseDistances::getClassType(void)
{
    
    static std::string rev_type = "Func_treePairwiseDistances";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_treePairwiseDistances::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<Tree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_treePairwiseDistances::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTreePairwiseDistances";
    
    return f_name;
}


const TypeSpec& Func_treePairwiseDistances::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

