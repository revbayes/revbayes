#include "NearestNeighborInterchangeProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TreeUtilities.h"
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
NearestNeighborInterchangeProposal::NearestNeighborInterchangeProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NearestNeighborInterchangeProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NearestNeighborInterchangeProposal* NearestNeighborInterchangeProposal::clone( void ) const
{
    
    return new NearestNeighborInterchangeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NearestNeighborInterchangeProposal::getProposalName( void ) const
{
    static std::string name = "NNI";
    
    return name;
}


double NearestNeighborInterchangeProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NearestNeighborInterchangeProposal::doProposal( void )
{
    
    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    if ( tau.getNumberOfTips() < 3)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick a random node which is not the root nor a direct descendant of the root
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
    
    // we need to work with the times
    double gparent_age  = grandparent.getAge();
    double parent_age   = parent.getAge();
    double uncles_age   = uncle->getAge();
    double my_age       = node->getAge();
    
    // now we store all necessary values
    storedChosenNode    = node;
    storedUncle         = uncle;
    
    storedAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*storedChosenNode, storedAges);
    TreeUtilities::getAges(*storedUncle, storedAges);

    // now exchange the two nodes
    grandparent.removeChild( uncle );
    parent.removeChild( node );
    grandparent.addChild( node );
    parent.addChild( uncle );
    node->setParent( &grandparent );
    uncle->setParent( &parent );
    
    // draw new ages and compute the hastings ratio at the same time
    // Note: the Hastings ratio needs to be there because one of the nodes might be a tip and hence not scaled!
    double lnHastingsratio = 0.0;
    
    // node rescale and MH
    double my_new_age;
    if ( node->isTip() )
    {
        my_new_age = my_age;
        lnHastingsratio = 0.0;
    }
    else
    {
        my_new_age = gparent_age * rng->uniform01();
        lnHastingsratio += log( gparent_age / parent_age);

        // rescale the subtrees
        double scaling_factor = my_new_age / my_age;
        TreeUtilities::rescaleSubtree(*node, scaling_factor );
        
        // compute the Hastings ratio
        size_t nNodes = node->getNumberOfNodesInSubtree(false);
        lnHastingsratio += (nNodes > 1 ? log( scaling_factor ) * (nNodes-1) : 0.0 );

    }

    // uncle rescale and MH
    double uncles_new_age;
    if ( uncle->isTip() )
    {
        uncles_new_age = uncles_age;
    }
    else
    {
        uncles_new_age = parent_age * rng->uniform01();
        lnHastingsratio += log( parent_age / gparent_age);
        
        double scaling_factor = uncles_new_age / uncles_age;
        TreeUtilities::rescaleSubtree(*uncle, scaling_factor );

        size_t nNodes = uncle->getNumberOfNodesInSubtree(false);
        lnHastingsratio += (nNodes > 1 ? log( scaling_factor ) * (nNodes-1) : 0.0 );
        
    }
    
    return lnHastingsratio;
    
}


/**
 *
 */
void NearestNeighborInterchangeProposal::prepareProposal( void )
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
void NearestNeighborInterchangeProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void NearestNeighborInterchangeProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        TopologyNode& parent = storedUncle->getParent();
        TopologyNode& grandparent = storedChosenNode->getParent();
    
        // now exchange the two nodes
        grandparent.removeChild( storedChosenNode );
        parent.removeChild( storedUncle );
        grandparent.addChild( storedUncle );
        parent.addChild( storedChosenNode );
        storedUncle->setParent( &grandparent );
        storedChosenNode->setParent( &parent );
    
        // rescale to old ages
        TreeUtilities::setAges(*storedChosenNode, storedAges);
        TreeUtilities::setAges(*storedUncle, storedAges);
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NearestNeighborInterchangeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NearestNeighborInterchangeProposal::setProposalTuningParameter(double tp)
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
void NearestNeighborInterchangeProposal::tune( double rate )
{
    
    // nothing to tune
    
}

