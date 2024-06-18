#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "DirichletTimeTreeDistribution.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;

DirichletTimeTreeDistribution::DirichletTimeTreeDistribution( const TypedDagNode<double> *r, const TypedDagNode<RbVector<double> > *a, const std::vector<Taxon> &n) : TypedDistribution<Tree>( new Tree() ),
    root_age( r ),
    alpha( a ),
    taxa( n )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( root_age );
    addParameter( alpha );

    num_taxa = taxa.size();
    
    simulateTree();
    
}


DirichletTimeTreeDistribution::~DirichletTimeTreeDistribution()
{
    
}


/**
 * Recursive call to attach ordered interior node times to the time tree psi. Call it initially with the
 * root of the tree.
 */
void DirichletTimeTreeDistribution::attachTimes(Tree* psi, std::vector<TopologyNode *> &nodes, size_t index, const std::vector<double> &interiorNodeTimes, double originTime )
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
        nodes.erase(nodes.begin()+long(node_index));
        
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
void DirichletTimeTreeDistribution::buildRandomBinaryHistory(std::vector<TopologyNode*> &tips)
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
        tips.erase(tips.begin()+long(index));
        
        // Add a left child
        TopologyNode* leftChild = new TopologyNode(0);
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // Add a right child
        TopologyNode* rightChild = new TopologyNode(0);
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // Recursive call to this function
        buildRandomBinaryHistory(tips);
    }
}


/* Clone function */
DirichletTimeTreeDistribution* DirichletTimeTreeDistribution::clone( void ) const
{
    
    return new DirichletTimeTreeDistribution( *this );
}


/* Compute probability */
double DirichletTimeTreeDistribution::computeLnProbability( void )
{
    
    // Variable declarations and initialization
    double ln_prob = 0.0;
    double age = root_age->getValue();
    
    // we need to check that the root age matches
    if ( fabs(age - value->getRoot().getAge() ) > 1E-7 )
    {
        return RbConstants::Double::neginf;
    }
    
    std::vector<double> internal_node_ages;
    
    // check nodes and get ages
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        // check that the ages are in correct chronological order
        // i.e., no child is older than its parent
        if ( the_node.isRoot() == false )
        {
            
            if ( (the_node.getAge() - the_node.getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( (the_node.getAge() - the_node.getParent().getAge()) > 1E-6 && the_node.isSampledAncestorTip() == true )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
        // check that the sampled ancestor nodes have a zero branch length
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
        
        if ( the_node.isTip() == false )
        {
            internal_node_ages.push_back( the_node.getAge() );
        }
        
    }
    
    // now sort the ages and the compute the differences
    std::sort(internal_node_ages.begin(), internal_node_ages.end());
    for (size_t i=internal_node_ages.size()-1; i>0; --i)
    {
        // take the difference between two consecutive ages
        internal_node_ages[i] -= internal_node_ages[i-1];
        // and then normalize
        internal_node_ages[i] /= age;
    }
    
    if ( alpha->getValue().size() != internal_node_ages.size() )
    {
        throw RbException("Missmatching size of alpha parameter (" + StringUtilities::toString(alpha->getValue().size()) + ") with number of between node age intervals (" +  StringUtilities::toString(internal_node_ages.size()) + ") in DirichletTimeTree distribution.");
    }
    
    // Take the Dirichlet draw into account
    ln_prob = RbStatistics::Dirichlet::lnPdf(alpha->getValue(), internal_node_ages);

    // Take the ordering effect into account
    ln_prob += RbMath::lnFactorial( int(num_taxa - 2) );

    return ln_prob;
}


void DirichletTimeTreeDistribution::redrawValue( void )
{
    simulateTree();
}


/** Simulate the tree conditioned on the time of origin */
void DirichletTimeTreeDistribution::simulateTree( void )
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
        nodes.erase(nodes.begin()+long(index) );
        
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

void DirichletTimeTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_age)
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void DirichletTimeTreeDistribution::keepSpecialization(const DagNode *affecter)
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
void DirichletTimeTreeDistribution::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->restoreAffected();
    }
    
}

/** Swap a parameter of the distribution */
void DirichletTimeTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == root_age)
    {
        root_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void DirichletTimeTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->touchAffected();
    }
    
}
