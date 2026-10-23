#include "NarrowExchangeProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
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
NarrowExchangeProposal::NarrowExchangeProposal( StochasticNode<Tree> *n, StochasticNode< RbVector<Tree> > *vec_n ) : Proposal(),
    variable( n ),
    vector_variable( vec_n )
{
    // tell the base class to add the node
    addNode( variable );
    addNode( vector_variable );

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NarrowExchangeProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NarrowExchangeProposal* NarrowExchangeProposal::clone( void ) const
{
    
    return new NarrowExchangeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NarrowExchangeProposal::getProposalName( void ) const
{
    static std::string name = "NarrowExchange";
    
    return name;
}


double NarrowExchangeProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new simplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NarrowExchangeProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // get the current tree either from the single variable or the vector of trees
    Tree *tmp = NULL;
    if ( variable != NULL )
    {
        tmp = &variable->getValue();
    }
    else
    {
        tree_index = floor(rng->uniform01() * vector_variable->getValue().size());
        tmp = &(vector_variable->getValue()[tree_index]);
    }
    Tree& tau = *tmp;
    
    if ( tau.getNumberOfTips() < 3)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick a random node which is not the root and neithor a direct descendant of the root
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
    if ( uncle == &parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    
    if (uncle->isSampledAncestorTip())
    {
        TopologyNode* brother = &parent.getChild(0);
        if (node == brother)
            brother = &parent.getChild(1);

        if (brother->isSampledAncestorTip())
        {
            // Swapping uncle and node, will mean that parent has two sampled-ancestor children.
            failed = true;
            return RbConstants::Double::neginf; // fail proposal
        }
    }

    // we need to work with the times
    double parent_age   = parent.getAge();
    double uncles_age   = uncle->getAge();

    // A sampled ancestor takes its destination parent's age. Allow it to move down as well
    // as up; otherwise an upward exchange has no reverse proposal.
    if (!uncle->isSampledAncestorTip() && uncles_age >= parent_age)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    // Work around setAge() refusing out-of-range assignments without throwing, and the tree
    // probability not enforcing fixed ages. Check before mutating so rejection leaves the tree intact.
    // These checks could be relaxed once assignment reliably completes and the tree probability
    // enforces all specified bounds. Equal bounds conservatively preserve the current fossil age.
    for (TopologyNode* fossil : {node, uncle})
    {
        if (!fossil->isSampledAncestorTip())
            continue;

        double new_age = fossil == node ? grandparent.getAge() : parent_age;
        const auto& interval = fossil->getTaxon().getAgeRange();
        if ((interval.zeroWidth() and new_age != fossil->getAge()) or
            (not interval.zeroWidth() and not interval.contains(new_age)))
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }

    failed = false;

    // now we store all necessary values
    storedChosenNode    = node;
    storedUncle         = uncle;

    // now exchange the two nodes
    grandparent.removeChild( uncle );
    parent.removeChild( node );
    grandparent.addChild( node );
    parent.addChild( uncle );
    node->setParent( &grandparent );
    uncle->setParent( &parent );
        
    if (node->isSampledAncestorTip())
        node->setAge( node->getParent().getAge());

    if (uncle->isSampledAncestorTip())
        uncle->setAge( uncle->getParent().getAge());

    // Selecting the exchanged uncle reverses the swap and restores sampled-ancestor ages.
    // Both selections have the same probability; internal ages and the eligible node count stay fixed.
    return 0.0;
}


/**
 *
 */
void NarrowExchangeProposal::prepareProposal( void )
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
void NarrowExchangeProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void NarrowExchangeProposal::undoProposal( void )
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
        
        if (storedChosenNode->isSampledAncestorTip())
            storedChosenNode->setAge( storedChosenNode->getParent().getAge());

        if (storedUncle->isSampledAncestorTip())
            storedUncle->setAge( storedUncle->getParent().getAge());
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NarrowExchangeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >(newN);
    }
    else if ( oldN == vector_variable )
    {
        vector_variable = static_cast<StochasticNode< RbVector<Tree> >* >( newN );
    }
    
}
