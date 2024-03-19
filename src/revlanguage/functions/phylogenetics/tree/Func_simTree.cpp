#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_simTree.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TopologyNode.h"

using namespace RevLanguage;

/** Default constructor */
Func_simTree::Func_simTree( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_simTree* Func_simTree::clone( void ) const
{
    
    return new Func_simTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_simTree::execute( void )
{
    
    int num_taxa             = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string& type = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    
    // the time tree object (topology + times)
    RevBayesCore::Tree *psi = new RevBayesCore::Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    RevBayesCore::TopologyNode* root = new RevBayesCore::TopologyNode();
    std::vector<RevBayesCore::TopologyNode* > nodes;
    nodes.push_back(root);
    
    if ( type == "balanced" )
    {
        simulateBalancedTree(num_taxa, nodes);
    }
    else if ( type == "caterpillar" )
    {
        simulateCaterpillarTree(num_taxa, root);
    }
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    // set the ages recursively
    setAges(psi, *root);
    
    return new RevVariable( new TimeTree( psi ) );
}


/** Get argument rules */
const ArgumentRules& Func_simTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "num_taxa", Natural::getClassTypeSpec(), "How many taxa this tree has.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "balanced" );
        optionsCondition.push_back( "caterpillar" );
//        optionsCondition.push_back( "random" );
        argumentRules.push_back( new OptionRule( "type"    , new RlString("balanced"), optionsCondition, "The type of the shape of the topology." ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_simTree::getClassType(void)
{
    
    static std::string rev_type = "Func_simTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_simTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_simTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_simTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_simTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    
    return return_typeSpec;
}


void Func_simTree::setAges(RevBayesCore::Tree *t, RevBayesCore::TopologyNode &n)
{
    
    if ( n.isTip() )
    {
        t->getNode( n.getIndex() ).setAge( 0.0 );
    }
    else
    {
        RevBayesCore::TopologyNode &left = n.getChild( 0 );
        setAges(t, left);
        
        RevBayesCore::TopologyNode &right = n.getChild( 1 );
        setAges(t, right);
        
        double a = t->getNode(left.getIndex()).getAge();
        double b = t->getNode(right.getIndex()).getAge();
        double max = (a > b ? a : b);
        
        t->getNode(n.getIndex()).setAge(max + 1.0);
    }
    
}


void Func_simTree::simulateBalancedTree( size_t n, std::vector<RevBayesCore::TopologyNode*> &nodes )
{
    
    // check if the number of taxa is divideable by two
    size_t half = n / 2;
    if ( (half+half) != n )
    {
        throw RbException("Bad number of taxa.");
    }
    
    std::vector<RevBayesCore::TopologyNode*> children;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        RevBayesCore::TopologyNode *parent = nodes[i];
        
        // add a left child
        RevBayesCore::TopologyNode* leftChild = new RevBayesCore::TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        children.push_back(leftChild);
        
        // add a right child
        RevBayesCore::TopologyNode* rightChild = new RevBayesCore::TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        children.push_back(rightChild);
        
    }
    
    if ( half == 1 )
    {
        // we are done with the recursion
        for (size_t i = 0; i < children.size(); ++i)
        {
            RevBayesCore::TopologyNode *node = children[i];
            std::string name = "Taxon_" + StringUtilities::to_string(i+1);
            node->setName(name);
        }
        
    }
    else
    {
        simulateBalancedTree(half, children);
    }
    
}




void Func_simTree::simulateCaterpillarTree( size_t n, RevBayesCore::TopologyNode* node )
{
    
    
    // add a left child
    RevBayesCore::TopologyNode* leftChild = new RevBayesCore::TopologyNode();
    node->addChild(leftChild);
    leftChild->setParent(node);
    
    // add a right child
    RevBayesCore::TopologyNode* rightChild = new RevBayesCore::TopologyNode();
    node->addChild(rightChild);
    rightChild->setParent(node);
    
    std::string name = "Taxon_" + StringUtilities::to_string(n);
    rightChild->setName(name);
    
    if ( n > 2 )
    {
        simulateCaterpillarTree(n-1, leftChild);
    }
    else
    {
        std::string name = "Taxon_" + StringUtilities::to_string(n-1);
        leftChild->setName(name);
    }
    
}

