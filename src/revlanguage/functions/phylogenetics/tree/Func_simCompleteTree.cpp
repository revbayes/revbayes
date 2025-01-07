#include <cmath>       /* log2 */
#include <cstddef>
#include <algorithm>
#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_simCompleteTree.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "RlTaxon.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Procedure.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** Default constructor */
Func_simCompleteTree::Func_simCompleteTree( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_simCompleteTree* Func_simCompleteTree::clone( void ) const
{
    
    return new Func_simCompleteTree( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_simCompleteTree::execute( void )
{
    
    size_t arg_index = 0;
    
    // taxon names
    const std::vector<RevBayesCore::Taxon> &taxa   = static_cast<const ModelVector<Taxon> &>( args[arg_index++].getVariable()->getRevObject() ).getDagNode()->getValue();

    double age  = static_cast<const RealPos &>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    double b    = static_cast<const RealPos &>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    double d    = static_cast<const RealPos &>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    // get the number of taxa
    size_t num_taxa = taxa.size();
    
    // the time tree object (topology + times)
    RevBayesCore::Tree *psi = new RevBayesCore::Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    std::vector<RevBayesCore::TopologyNode* > nodes = std::vector<RevBayesCore::TopologyNode* >(num_taxa, NULL);
    for (size_t i=0; i<num_taxa; ++i)
    {
        
        // create the i-th taxon
        RevBayesCore::TopologyNode* node = new RevBayesCore::TopologyNode( taxa[i], i );
        
        // set the age of this tip node
        node->setAge( taxa[i].getAge() );
        
        // add the new node to the list
        nodes[i] = node;
    }
    // Get the rng
    RevBayesCore::RandomNumberGenerator* rng = RevBayesCore::GLOBAL_RNG;
    

    std::vector<double> times = std::vector<double>(num_taxa-1, 0.0);
    for (size_t i=0; i<(num_taxa-2); ++i)
    {
        
        double rho = 1.0;
    
        // get a random draw
        double u = rng->uniform01();
        
        // compute the time for this draw
        double t = 0.0;
        if ( b > d )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(rho*b+(b*(1-rho)-d)*exp((d-b)*age) ) ) ) - (b*(1-rho)-d) ) / (rho * b) ) + (d-b)*age )  /  (d-b);
        }
        else
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(rho*b*exp((b-d)*age)+(b*(1-rho)-d) ) ) ) - (b*(1-rho)-d) ) / (rho * b) ) + (d-b)*age )  /  (d-b);
        }
    
        times[i] = age - t;
    }
    times[num_taxa-2] = age;
    
    // finally sort the times
    std::sort(times.begin(), times.end());
    
    buildRandomBinaryTree(nodes, times);
    
    RevBayesCore::TopologyNode *root = nodes[0];
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    return new RevVariable( new TimeTree( psi ) );
}


void Func_simCompleteTree::buildRandomBinaryTree(std::vector<RevBayesCore::TopologyNode*> &nodes, std::vector<double> &ages)
{
    
    if ( nodes.size() > 1 )
    {
        // Get the rng
        RevBayesCore::RandomNumberGenerator* rng = RevBayesCore::GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*nodes.size()) );
        
        // get the node from the list
        RevBayesCore::TopologyNode* left_child = nodes[index];
        
        // remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+long(index));
        
        
        
        // repeat
        index = static_cast<size_t>( floor(rng->uniform01()*nodes.size()) );
        
        // get the node from the list
        RevBayesCore::TopologyNode* right_child = nodes[index];
        
        // remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+long(index));
        
        
        RevBayesCore::TopologyNode* parent = new RevBayesCore::TopologyNode();
        parent->setAge( ages[0] );
        ages.erase( ages.begin() );
        
        // add a left child
        parent->addChild(left_child);
        left_child->setParent(parent);
        
        // add a right child
        parent->addChild(right_child);
        right_child->setParent(parent);
        
        
        nodes.push_back( parent );
        
        // recursive call to this function
        buildRandomBinaryTree( nodes, ages );
    }
    
}


/** Get argument rules */
const ArgumentRules& Func_simCompleteTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "taxa"    , ModelVector<Taxon>::getClassTypeSpec(), "The taxon names used for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "rootAge" , RealPos::getClassTypeSpec(), "The root age.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "lambda"  , RealPos::getClassTypeSpec(), "The speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mu"      , RealPos::getClassTypeSpec(), "The extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_simCompleteTree::getClassType(void)
{
    
    static std::string rev_type = "Func_simCompleteTree";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_simCompleteTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_simCompleteTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simCompleteTree";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_simCompleteTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_simCompleteTree::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = TimeTree::getClassTypeSpec();
    
    return return_typeSpec;
}
