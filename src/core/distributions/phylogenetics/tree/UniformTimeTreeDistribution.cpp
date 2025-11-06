#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "UniformTimeTreeDistribution.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;

UniformTimeTreeDistribution::UniformTimeTreeDistribution( const TypedDagNode<double> *a, const std::vector<Taxon> &n) : TypedDistribution<Tree>( new Tree() ),
    root_age( a ),
    taxa( n )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( root_age );
    
    num_taxa = taxa.size();
    
    simulateTree();
    
}


UniformTimeTreeDistribution::~UniformTimeTreeDistribution()
{
    
}


/**
 * Recursive call to attach ordered interior node times to the time tree psi. Call it initially with the
 * root of the tree.
 */
void UniformTimeTreeDistribution::attachTimes(Tree* psi, std::vector<TopologyNode *> &nodes, size_t index, const std::vector<double> &interiorNodeTimes, double originTime )
{
    
    if (index < num_taxa-1)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // Randomly draw one node from the list of nodes
        size_t node_index = static_cast<size_t>( floor(rng->uniform01()*nodes.size()) );
        
        // Get the node from the list
        TopologyNode* parent = nodes.at(node_index);
        psi->getNode( parent->getIndex() ).setAge( originTime - interiorNodeTimes[index] );
        
        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+std::int64_t(node_index));
        
        // Add the left child if an interior node
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            nodes.push_back(leftChild);
        }

        // Add the right child if an interior node
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            nodes.push_back(rightChild);
        }
        
        // Recursive call to this function
        attachTimes(psi, nodes, index+1, interiorNodeTimes, originTime);
    }
    
}


/** Build random binary tree to size num_taxa. The result is a draw from the uniform distribution on histories. */
void UniformTimeTreeDistribution::buildRandomBinaryHistory(std::vector<TopologyNode*> &tips)
{
    
    if (tips.size() < num_taxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;

        // Randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );

        // Get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // Remove the randomly drawn node from the list
        tips.erase(tips.begin()+std::int64_t(index));
        
        // Add a left child
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // Add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // Recursive call to this function
        buildRandomBinaryHistory(tips);
    }
}


/* Clone function */
UniformTimeTreeDistribution* UniformTimeTreeDistribution::clone( void ) const
{
    
    return new UniformTimeTreeDistribution( *this );
}


/* Compute probability */
double UniformTimeTreeDistribution::computeLnProbability( void )
{
    
    // Variable declarations and initialization
    double lnProb = 0.0;
    double age = root_age->getValue();
    
    // we need to check that the root age matches
    if ( fabs(age - value->getRoot().getAge() ) > 1E-7 )
    {
        return RbConstants::Double::neginf;
    }
    
    
    // check that the ages are in correct chronological order
    // i.e., no child is older than its parent
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        if ( the_node.isRoot() == false )
        {
            
            if ( (the_node.getAge() - (*it)->getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( (the_node.getAge() - (*it)->getParent().getAge()) > 1E-6 && the_node.isSampledAncestorTip() == true )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    // check that the sampled ancestor nodes have a zero branch length
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        if ( the_node.isSampledAncestorTip() == true )
        {
            
            if ( the_node.isFossil() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( the_node.getBranchLength() > 1E-6 )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    // Take the uniform draws into account
    lnProb = (num_taxa - 2) * log( 1.0 / age );

    // Take the ordering effect into account
    lnProb += RbMath::lnFactorial( int(num_taxa - 2) );

    return lnProb;
}


void UniformTimeTreeDistribution::redrawValue( void )
{
    simulateTree();
}


/** Simulate the tree conditioned on the time of origin */
void UniformTimeTreeDistribution::simulateTree( void )
{
    
    delete value;
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Create the time tree object (topology + times)
    Tree *psi = new Tree();

    // Root the topology by setting the appropriate flag
    psi->setRooted( true );
    
    // Create the root node and a vector of nodes
    TopologyNode* root = new TopologyNode();
    std::vector<TopologyNode* > nodes;
    nodes.push_back(root);
    
    // Draw a random tree history
    buildRandomBinaryHistory(nodes);
    
    // Set the tip names
    for (size_t i=0; i<num_taxa; i++)
    {
        size_t index = size_t( floor(rng->uniform01() * nodes.size()) );
        
        // Get the node from the list
        TopologyNode* node = nodes.at(index);
        
        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+std::int64_t(index) );
        
        // Set taxon
        node->setTaxon( taxa[i] );
    }
    
    // Initialize the topology by setting the root
    psi->setRoot(root, true);
    
    // Now simulate the speciation times counted from originTime
    std::vector<double> intNodeTimes;
    double              t_or = root_age->getValue();
    intNodeTimes.push_back( 0.0 );  // For time of mrca
    for ( size_t i=0; i<num_taxa-2; i++ )
    {
        double t = rng->uniform01() * t_or;
        intNodeTimes.push_back( t );
    }
    
    // Sort the speciation times from 0.0 (root node) to the largest value
    std::sort( intNodeTimes.begin(), intNodeTimes.end() );

    // Attach times
    nodes.clear();
    nodes.push_back( root );
    attachTimes(psi, nodes, 0, intNodeTimes, t_or);
    for (size_t i = 0; i < num_taxa; ++i)
    {
        TopologyNode& node = psi->getTipNodeWithName(taxa[i].getName());
        node.setAge( 0.0 );
    }
    
    // Finally store the new value
    value = psi;
    
}

void UniformTimeTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_age)
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void UniformTimeTreeDistribution::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age )
    {
        dag_node->keepAffected();
    }
    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UniformTimeTreeDistribution::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->restoreAffected();
    }
    
}

/** Swap a parameter of the distribution */
void UniformTimeTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == root_age)
    {
        root_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UniformTimeTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->touchAffected();
    }
    
}
