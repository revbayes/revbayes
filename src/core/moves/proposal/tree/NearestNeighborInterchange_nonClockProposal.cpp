#include "NearestNeighborInterchange_nonClockProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NearestNeighborInterchange_nonClockProposal::NearestNeighborInterchange_nonClockProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NearestNeighborInterchange_nonClockProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NearestNeighborInterchange_nonClockProposal* NearestNeighborInterchange_nonClockProposal::clone( void ) const
{
    
    return new NearestNeighborInterchange_nonClockProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NearestNeighborInterchange_nonClockProposal::getProposalName( void ) const
{
    static std::string name = "NNI";
    
    return name;
}


double NearestNeighborInterchange_nonClockProposal::getProposalTuningParameter( void ) const
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
double NearestNeighborInterchange_nonClockProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    TopologyNode& parent = node->getParent();
    TopologyNode& grandparent = parent.getParent();
    TopologyNode* uncle = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == &parent ) {
        uncle = &grandparent.getChild( 1 );
    }
    
    // now we store all necessary values
    storedChoosenNode   = node;
    storedUncle         = uncle;
    
    // now exchange the two nodes
    grandparent.removeChild( uncle );
    parent.removeChild( node );
    grandparent.addChild( node );
    parent.addChild( uncle );
    node->setParent( &grandparent );
    uncle->setParent( &parent );
    
    return 0.0;
    
}


/**
 *
 */
void NearestNeighborInterchange_nonClockProposal::prepareProposal( void )
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
void NearestNeighborInterchange_nonClockProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NearestNeighborInterchange_nonClockProposal::undoProposal( void )
{
    // undo the proposal
    TopologyNode& parent = storedUncle->getParent();
    TopologyNode& grandparent = storedChoosenNode->getParent();
    
    // now exchange the two nodes
    grandparent.removeChild( storedChoosenNode );
    parent.removeChild( storedUncle );
    grandparent.addChild( storedUncle );
    parent.addChild( storedChoosenNode );
    storedUncle->setParent( &grandparent );
    storedChoosenNode->setParent( &parent );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NearestNeighborInterchange_nonClockProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NearestNeighborInterchange_nonClockProposal::setProposalTuningParameter(double tp)
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
void NearestNeighborInterchange_nonClockProposal::tune( double rate )
{
    
    // nothing to tune
    
}

