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
    
    // we need to work with the times
    double parent_age   = parent.getAge();
    double uncles_age   = uncle->getAge();
    
    if ( uncles_age < parent_age )
    {
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
        
        return 0.0;
    }
    else
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
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


void NarrowExchangeProposal::setProposalTuningParameter(double tp)
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
void NarrowExchangeProposal::tune( double rate )
{
    
    // nothing to tune
    
}

