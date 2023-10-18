#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTree.h"
#include "RlTimeTree.h"
#include "RlUtils.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** Default constructor */
BranchLengthTree::BranchLengthTree(void) : Tree()
{
    
    initMethods();

}

/** Construct from bool */
BranchLengthTree::BranchLengthTree(RevBayesCore::Tree *t) : Tree( t )
{
    
    initMethods();

}

/** Construct from bool */
BranchLengthTree::BranchLengthTree(const RevBayesCore::Tree &t) : Tree( new RevBayesCore::Tree( t ) )
{
    
    initMethods();

}

/** Construct from core type */
BranchLengthTree::BranchLengthTree(RevBayesCore::TypedDagNode<RevBayesCore::Tree> *n) : Tree( n )
{
    
    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
BranchLengthTree* BranchLengthTree::clone(void) const
{

	return new BranchLengthTree(*this);
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> BranchLengthTree::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    if (name == "makeUltrametric")
    {

        found = true;

        RevBayesCore::Tree *tree = dag_node->getValue().clone();
        RevBayesCore::TreeUtilities::makeUltrametric(*tree);
        
        return new RevVariable( new TimeTree( tree ) );
    }
    
    return Tree::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& BranchLengthTree::getClassType(void)
{

    static std::string rev_type = "BranchLengthTree";

	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& BranchLengthTree::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Tree::getClassTypeSpec() ) );

	return rev_type_spec;
}


/** Get type spec */
const TypeSpec& BranchLengthTree::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
/**
 * Initialize the member methods.
 */
void BranchLengthTree::initMethods( void )
{
    ArgumentRules* rerootArgRules = new ArgumentRules();
    rerootArgRules->push_back( new ArgumentRule("leaf", RlString::getClassTypeSpec(), "The outgroup leaf.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction( new MemberProcedure( "reroot", RlUtils::Void,  rerootArgRules       ) );

    ArgumentRules* makeUltraArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "makeUltrametric", RlUtils::Void, makeUltraArgRules ) );
}
