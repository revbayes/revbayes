#include "BiasedNearestNeighborInterchangeProposal.h"

#include <stddef.h>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BiasedNearestNeighborInterchangeProposal::BiasedNearestNeighborInterchangeProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    attempts( 0 )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BiasedNearestNeighborInterchangeProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BiasedNearestNeighborInterchangeProposal* BiasedNearestNeighborInterchangeProposal::clone( void ) const
{
    
    return new BiasedNearestNeighborInterchangeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BiasedNearestNeighborInterchangeProposal::getProposalName( void ) const
{
    static std::string name = "NNI";
    
    return name;
}


double BiasedNearestNeighborInterchangeProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A SPR proposal.
 *
 * \return The hastings ratio.
 */
double BiasedNearestNeighborInterchangeProposal::doProposal( void )
{
    
    // add the current tree
    if ( attempts >= attempts_before_learning && attempts < attempts_to_learning )
    {
        storeTree();
    }
    
    // increase our attempts counter
    ++attempts;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random internal branch on which we are going to perform the NNI
    // the branch is represented by it's tipward node, so we can pick any internal node
    // which is not the root node
    TopologyNode* node;
    
    double forward_prob   = 0.0;
    double backwards_prob = 0.0;
    
    // check if we should use the biased proposal
    if ( attempts > attempts_before_using )
    {
        std::map<RbBitSet, TopologyNode*> current_clades = tau.getBitsetToNodeMap();
        size_t num_clades = current_clades.size();
        
        std::vector<double>        this_clade_frequencies = std::vector<double>(num_clades,0);
        std::vector<TopologyNode*> this_clade_nodes       = std::vector<TopologyNode*>(num_clades,NULL);
        double total_probability = 0.0;
        size_t xxx = 0;
        for ( std::map<RbBitSet, TopologyNode*>::const_iterator it=current_clades.begin(); it!=current_clades.end(); ++it)
        {
            double this_freq = clade_frequencies[it->first];
            this_clade_frequencies[xxx] = 1.0 - double(this_freq) / attempts;
            this_clade_nodes[xxx] = it->second;
            
            total_probability += 1.0 - double(this_freq) / attempts;
            
            ++xxx;
        }
        
        double u = rng->uniform01() * total_probability;
        xxx = 0;
        while ( u > this_clade_frequencies[xxx] )
        {
            u -= this_clade_frequencies[xxx];
            ++xxx;
        }
        node = this_clade_nodes[xxx];
        forward_prob = log( this_clade_frequencies[xxx] / total_probability );
    }
    else
    {
        do {
            double u = rng->uniform01();
            size_t index = size_t( std::floor((tau.getNumberOfNodes()-tau.getNumberOfTips()) * u) ) + tau.getNumberOfTips();
            node = &tau.getNode(index);
        } while ( node->isRoot() || node->isTip() );
    }
    
    // some important nodes that we need several times
    TopologyNode& parent = node->getParent();
    
    // some flags to tell us how to perform the tree modifications
    picked_uncle = true;
    picked_root_branch = parent.isRoot();
    
    // we first need to decide if we use the left or right child of this node
    TopologyNode* node_A = NULL;
    if ( rng->uniform01() < 0.5 )
    {
        node_A = &node->getChild(0);
    }
    else
    {
        node_A = &node->getChild(1);
    }
    
    // now we randomly pick the second node for the swap
    // here we need to check if this branch "begins" at the root
    TopologyNode* node_B = NULL;
    if ( picked_root_branch == false )
    {
        // it does not start at the root
        // so we can either pick the sibling, or the grand parent
        if ( rng->uniform01() < 0.5 )
        {
            // we picked the sibling of the chose node
            node_B = &parent.getChild(0);
            if ( node_B == node )
            {
                node_B = &parent.getChild(1);
            }
            
            picked_uncle = true;
        }
        else
        {
            // we picked the grant parent
            node_B = &parent.getParent();
            
            picked_uncle = false;
        }
    }
    else
    {
        // the branch begins at the root
        
        // first start by checking if the chosen node is child 0
        size_t child_offset = (node == &parent.getChild(0) ? 1 : 0);
        
        // now we can pick the second node, which can be either of the two remaining children
        if ( rng->uniform01() < 0.5 )
        {
            // we picked the first remaining child
            node_B = &parent.getChild(child_offset);
        }
        else
        {
            // we picked the second remaining child
            node_B = &parent.getChild(child_offset+1);
            if ( node_B == node )
            {
                node_B = &parent.getChild(child_offset+2);
            }
        }
        
        picked_uncle = true;
        
    }

    // now we store all necessary values
    stored_node_A = node_A;
    stored_node_B = node_B;

    if ( picked_root_branch == false && picked_uncle == true )
    {
        
        // now exchange the two nodes
        parent.removeChild( node_B );
        node->removeChild( node_A );
        parent.addChild( node_A );
        node->addChild( node_B );
        node_A->setParent( &parent );
        node_B->setParent( node );
    }
    else if ( picked_root_branch == false && picked_uncle == false )
    {
        node_B->removeChild( &parent );
        parent.removeChild( node );
        node->removeChild( node_A );
        node_B->addChild( node );
        node->addChild( &parent );
        parent.addChild( node_A );
        node->setParent( node_B );
        parent.setParent( node );
        node_A->setParent( &parent );
        
    }
    else if ( picked_root_branch == true  )
    {
        
        // now exchange the two nodes
        parent.removeChild( node_B );
        node->removeChild( node_A );
        parent.addChild( node_A );
        node->addChild( node_B );
        node_A->setParent( &parent );
        node_B->setParent( node );
    }
    
    
    if ( attempts > attempts_before_using )
    {
        std::map<RbBitSet, TopologyNode*> current_clades = tau.getBitsetToNodeMap();
        size_t num_clades = current_clades.size();
        
        double total_probability = 0.0;
        double chosen_probability = 0.0;
        size_t xxx = 0;
        for ( std::map<RbBitSet, TopologyNode*>::const_iterator it=current_clades.begin(); it!=current_clades.end(); ++it)
        {
            double this_freq = clade_frequencies[it->first];
            this_clade_frequencies[xxx] = 1.0 - double(this_freq) / attempts;
            if ( it->second == store);
            
            total_probability += 1.0 - double(this_freq) / attempts;
            
            ++xxx;
        }
        
        backwards_prob = log( chosen_probability / total_probability );
    }
    
    return backwards_prob - forward_prob;
}


/**
 *
 */
void BiasedNearestNeighborInterchangeProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void BiasedNearestNeighborInterchangeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}



/**
 * Store the current tree to our list of observed bipartition frequencies.
 */
void BiasedNearestNeighborInterchangeProposal::storeTree( void )
{
    
    // get the current tree
    const Tree& my_tree = tree->getValue();
        
    // get the number of tips for this tree to create bitsets of the correct size
    size_t num_tips = my_tree.getNumberOfTips();
    
    // create our vector of clades (represented as bitsets)
    std::vector<RbBitSet> all_clades;
    
    // get the root node to start the clade collection
    const TopologyNode& root = my_tree.getRoot();

    // get the children of the root node
    const std::vector<TopologyNode* > root_children = root.getChildren();
    
    // now extract all the clades as bitsets for all nodes
    for ( std::vector<TopologyNode* >::const_iterator i=root_children.begin(); i!=root_children.end(); ++i )
    {
        (*i)->getAllClades(all_clades, num_tips, true);
    }
    
    for ( std::vector<RbBitSet>::const_iterator i=all_clades.begin(); i!=all_clades.end(); ++i )
    {
        const RbBitSet& this_clade = *i;
        RbBitSet this_clade_complement = this_clade;
        this_clade_complement.flip();
        
        std::map<RbBitSet,size_t>::const_iterator this_clade_freq = clade_frequencies.find( this_clade );
        
        // check if we have seen this clade before
        if ( this_clade_freq == clade_frequencies.end() )
        {
            clade_frequencies.insert( std::make_pair(this_clade, 1) );
            clade_frequencies.insert( std::make_pair(this_clade_complement, 1) );
        }
        else
        {
            ++clade_frequencies[this_clade];
            ++clade_frequencies[this_clade_complement];
        }
    }
    
    
}



/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void BiasedNearestNeighborInterchangeProposal::undoProposal( void )
{
    // undo the proposal
    TopologyNode* parent = &stored_node_A->getParent();
    TopologyNode* node   = ( (picked_root_branch == false && picked_uncle == false) ? &parent->getParent() : &stored_node_B->getParent());
    
    TopologyNode* node_A = stored_node_A;
    TopologyNode* node_B = stored_node_B;
    
    // now exchange the two nodes
    if ( picked_root_branch == false && picked_uncle == true )
    {
        
        // now exchange the two nodes
        parent->removeChild( node_A );
        node->removeChild( node_B );
        parent->addChild( node_B );
        node->addChild( node_A );
        node_A->setParent( node );
        node_B->setParent( parent );
    }
    else if ( picked_root_branch == false && picked_uncle == false )
    {
        node = &parent->getParent();
            
        node_B->removeChild( node );
        parent->removeChild( node_A );
        node->removeChild( parent );
        node_B->addChild( parent );
        node->addChild( node_A );
        parent->addChild( node );
        node->setParent( parent );
        parent->setParent( node_B );
        node_A->setParent( node );
        
    }
    else if ( picked_root_branch == true  )
    {
        
        // now exchange the two nodes
        parent->removeChild( node_A );
        node->removeChild( node_B );
        parent->addChild( node_B );
        node->addChild( node_A );
        node_A->setParent( node );
        node_B->setParent( parent );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BiasedNearestNeighborInterchangeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void BiasedNearestNeighborInterchangeProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void BiasedNearestNeighborInterchangeProposal::tune( double rate )
{
    
    // nothing to tune
    
}

