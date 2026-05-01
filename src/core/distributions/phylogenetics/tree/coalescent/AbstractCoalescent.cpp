#include "AbstractCoalescent.h"

#include <cmath>
#include <iosfwd>

#include "Clade.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


/**
 * Default Constructor.
 *
 * The constructor connects the parameters of the coalescent process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 *
 * \param[in] tn        A vector of taxon names used during initialization.
 * \param[in] c         a vector of type Clade for clade constraints.
 */
AbstractCoalescent::AbstractCoalescent(const std::vector<Taxon> &tn, const std::vector<Clade> &c) : TypedDistribution<Tree>( new Tree() ),
    constraints( c ),
    num_taxa( tn.size() ),
    taxa( tn )
{
    
    // the combinatorial factor for the probability of a labelled history is
    // 2^{n-1} / ( n! * (n-1)! )
    // but this probability cancels with sum( choose(k,2) ) and thus we ommit the tree probability
    // You need to make sure that we cancel this out in the probability of each coalescent time!!
    
    logTreeTopologyProb = 0.0;

}



/**
 * Randomly attach the node times to a tree topology.
 * This function works by randomly picking a node from the set of tips,
 * setting its time to times[index], increment the index,
 * adding the two children (if they are not actual tips) to the set of tips,
 * and recursively calling this function again.
 *
 * \param[in]     psi        The tree topology (needed to call setAge).
 * \param[in]     tips       The vector of TopologyNode for the tips
 * \param[in]     index      The index of the node to be aged. Used for recursion
 * \param[in]     ages       A vector of node ages
 */
void AbstractCoalescent::attachAges(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &ages)
{
    
    if ( index < ages.size() )
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t tip_index = size_t( floor(rng->uniform01()*tips.size()) );
        
        // get the node from the list
        TopologyNode* parent = tips.at(tip_index);
        psi->getNode( parent->getIndex() ).setAge( ages[ages.size()-index-1] );
        
        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+ long(tip_index) );
        
        // add a left child
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            tips.push_back(leftChild);
        }
        
        // add a right child
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            tips.push_back(rightChild);
        }
        
        // recursive call to this function
        attachAges(psi, tips, index+1, ages);
    }
    
    // no return value
}

/*
 * Builds a Random Binary Tree by starting with tips and coalescening them together randomly
 *
 * @param tips a vector of TopologyNode for the tips
 *
 */
void AbstractCoalescent::buildRandomBinaryTree(std::vector<TopologyNode*> &tips)
{
    
    if (tips.size() < num_taxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
        
        // get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+long(index));
        
        // add a left child
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // recursive call to this function
        buildRandomBinaryTree(tips);
    }
    
}

/**
 * Randomly build a tree by picking two tips and coalescing, while also attaching the times to a tree topology.
 * This function works by randomly picking two node from the set of tips, adding their parent,
 * setting the parent's time is set to ages[ageIndex], incrementing ageIndex and index (for node numbering),
 * and recursively calling this function again.
 *
 * \param[in]     psi        The tree topology (needed to call setAge).
 * \param[in]     active     The vector of nodes eligible to be coalesced
 * \param[in]     ages       The ages of the coalescent events
 */
void AbstractCoalescent::buildHeterochronousRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &active, const std::vector<double> &ages)
{
    for (size_t i = 0; i < (num_taxa - 1); ++i)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        bool valid = false;
        size_t redraw_attempts = 0;
        std::vector<TopologyNode* > unusableNodes;
        do
        {
            
            // randomly draw one child (arbitrarily called left) node from the list of active nodes
            size_t left = static_cast<size_t>( floor(rng->uniform01()*active.size()) );
            TopologyNode* leftChild = active.at(left);
            
            // remove the randomly drawn node from the list
            active.erase(active.begin()+long(left));
            
            // randomly draw one child (arbitrarily called left) node from the list of active nodes
            size_t right = static_cast<size_t>( floor(rng->uniform01()*active.size()) );
            TopologyNode* rightChild = active.at(right);
            
            // remove the randomly drawn node from the list
            active.erase(active.begin()+long(right));
            
            // check that we aren't trying to coalesce nodes before they exist
            if ( (leftChild->getAge() < ages[i]) && (rightChild->getAge() < ages[i]))
            {
                // add the parent
                TopologyNode* parent = new TopologyNode(i + num_taxa);
                parent->addChild(leftChild);
                parent->addChild(rightChild);
                leftChild->setParent(parent);
                rightChild->setParent(parent);
                parent->setAge( ages[i] );
                active.push_back(parent);
                
                // we coalesced, we can move on now
                valid = true;
            }
            else
            {
                ++redraw_attempts;
                // abort
                // if the node is coalesceable for this age, we put it back in active
                // otherwise, we collect it in unusableNodes so we can't try to re-draw it
                
                // handle left child
                if (leftChild->getAge() < ages[i])
                {
                    // this is node is currently coalesceable
                    active.push_back(leftChild);
                }
                else
                {
                    // this node cannot be coalesced right now
                    unusableNodes.push_back(leftChild);
                }
                // handle left child
                if (rightChild->getAge() < ages[i])
                {
                    // this is node is currently coalesceable
                    active.push_back(rightChild);
                }
                else
                {
                    // this node cannot be coalesced right now
                    unusableNodes.push_back(rightChild);
                }
            }
        } while ( !valid );
        
        // If we had to temporarily discard nodes, now we put them back
        if (redraw_attempts > 0)
        {
            for (size_t index = 0; index < unusableNodes.size(); ++index)
            {
                TopologyNode* thisNode = unusableNodes.at(index);
                active.push_back(thisNode);
            }
        }
    }

}

/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * @return double for the log-transformed probability
 *
 */
double AbstractCoalescent::computeLnProbability( void )
{
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    
    // first check if the current tree matches the clade constraints
    if ( !matchesConstraints() )
    {
        return RbConstants::Double::neginf;
    }
    
    // check that all children are younger than there parents
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        if ( !(*it)->isRoot() )
        {
            if ( (*it)->getAge() >= (*it)->getParent().getAge() )
            {
                return RbConstants::Double::neginf;
            }
        }
    }
    
    // multiply the probability of a descendant of the initial species
    lnProbTimes += computeLnProbabilityTimes();
    
    return lnProbTimes + logTreeTopologyProb;
    
}



/**
 * We check here if all the constraints are satisfied.
 * These are hard constraints, that is, the clades must be monophyletic.
 *
 * \return     True if the constraints are matched, false otherwise.
 */
bool AbstractCoalescent::matchesConstraints( void )
{
    
    const TopologyNode &root = value->getRoot();
    
    for (std::vector<Clade>::iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        if ( !root.containsClade( *it, true ) )
        {
            return false;
        }
    }
    
    return true;
}


/**
 * Redraw the current value. We delegate this to the simulate method.
 */
void AbstractCoalescent::redrawValue( void )
{

    simulateTree();

}


/**
 *
 */
void AbstractCoalescent::simulateTree( void )
{
    bool homochronous = true;
    for (size_t i=0; i<num_taxa; ++i)
    {
        if ( taxa[i].getAge() > 0.0 )
        {
          homochronous = false;
          break;
        }
    }
    if ( homochronous )
    {
        simulateHomochronousTree();
    }
    else
    {
        simulateHeterochronousTree();
    }
    
}

/**
 *
 */
void AbstractCoalescent::simulateHeterochronousTree( void )
{
    
    // the time tree object (topology + times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    // make a vector of tip nodes
    std::vector<TopologyNode* > nodes;

    // set tip names
    for (size_t i=0; i<num_taxa; ++i)
    {
        // get the node from the list
        TopologyNode* node = new TopologyNode(i);
        
        // set name and age
        const std::string& name = taxa[i].getName();
        node->setName(name);
        node->setSpeciesName(taxa[i].getSpeciesName());
        node->setAge(taxa[i].getAge());
        // add to tips
        nodes.push_back(node);
    }
    
    // get times for simulation
    std::vector<double> ages = simulateCoalescentAges(num_taxa-1);
    
    // recursively build the tree
    buildHeterochronousRandomBinaryTree(psi, nodes, ages);

    // initialize the topology by setting the root
    // the root is the only node left in nodes
    TopologyNode* root = nodes[0]; // Only node left after coalescing all is the root
    psi->setRoot(root, true);

    // finally store the new value
    delete value;
    value = psi;
    
}

/** Simulates a coalescent tree
 *
 *
 *@throw RbException if the vector of taxan has nodes with different ages. Taxa cannot be heterochronous.
 */
void AbstractCoalescent::simulateHomochronousTree( void )
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // the time tree object (topology + times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    TopologyNode* root = new TopologyNode();
    std::vector<TopologyNode* > nodes;
    nodes.push_back(root);
    
    // recursively build the tree
    buildRandomBinaryTree(nodes);
    
    // set tip names
    for (size_t i=0; i<num_taxa; ++i)
    {
        size_t index = size_t( floor(rng->uniform01() * nodes.size()) );
        
        // get the node from the list
        TopologyNode* node = nodes.at(index);
        
        // remove the randomly drawn node from the list
        nodes.erase( nodes.begin()+long(index) );
        
        // set name
        const std::string& name = taxa[i].getName();
        node->setName(name);
        node->setSpeciesName(taxa[i].getSpeciesName());

    }
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    nodes.clear();
    
    if ( 1 < num_taxa)
    {
        // draw a time for each coalescent event
        std::vector<double> ages = simulateCoalescentAges(num_taxa-1);
        
        // add a left child
        TopologyNode* leftChild = &root->getChild(0);
        if ( !leftChild->isTip() )
        {
            nodes.push_back(leftChild);
        }
            
        // add a right child
        TopologyNode* rightChild = &root->getChild(1);
        if ( !rightChild->isTip() )
        {
            nodes.push_back(rightChild);
        }
        
        attachAges(psi, nodes, 1, ages);
        
        psi->getNode( root->getIndex() ).setAge( ages[ages.size()-1]);
        
    }
    
    // \todo Why are we doing this? (Sebastian)
    for (size_t i = 0; i < num_taxa; ++i)
    {
        TopologyNode& node = psi->getTipNode(i);
        psi->getNode( node.getIndex() ).setAge( 0.0 );
    }
    
    // finally store the new value
    delete value;
    value = psi;
    
}

