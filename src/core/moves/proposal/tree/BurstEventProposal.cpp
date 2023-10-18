#include <cstddef>
#include <cmath>
#include <iostream>

#include "BirthDeathBurstProcess.h"
#include "BurstEventProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BurstEventProposal::BurstEventProposal( StochasticNode<Tree> *n ) : Proposal(),
    psi( n )
{

    // tell the base class to add the node
    addNode( psi );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BurstEventProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BurstEventProposal* BurstEventProposal::clone( void ) const
{
    
    return new BurstEventProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BurstEventProposal::getProposalName( void ) const
{
    static std::string name = "BurstEvent";
    
    return name;
}


double BurstEventProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double BurstEventProposal::doProposal( void )
{
    
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree &tree = psi->getValue();
    BirthDeathBurstProcess *dist = static_cast<BirthDeathBurstProcess *>( &psi->getDistribution() );
    
    double burst_age = dist->getBurstTime();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    size_t node_index;
    do {
        double u = rng->uniform01();
        node_index = size_t( std::floor(tree.getNumberOfNodes() * u) );
        node = &tree.getNode(node_index);
    } while ( node->isRoot() || node->isTip() );
    
    TopologyNode& parent = node->getParent();
    TopologyNode& child = node->getChild(0);
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = child.getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child = node->getChild(1);
        child_Age = child.getAge();
    }

    if ( parent_age < burst_age || child_Age > burst_age ||
            dist->isBurstSpeciation(parent.getIndex()) ||
            dist->isBurstSpeciation(child.getIndex()) )
    {
        failed = true;
        // TODO: Should this be returning -inf (i.e., always reject proposal)?
        return 0;
    }
    
    // store variable
    stored_node = node;
    stored_age = my_age;
    
    double hr = 0;
    if ( dist->isBurstSpeciation(node_index) == true )
    {
        
        // draw new ages and compute the hastings ratio at the same time
        double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;
        
        // set the age
        stored_node->setAge( my_new_age );
        
        dist->setBurstSpeciation(node_index, false);
        
        hr += log(parent_age - child_Age);
    }
    else
    {
        // set the age
        stored_node->setAge( burst_age );
        
        dist->setBurstSpeciation(node_index, true);
        
        hr += - log(parent_age - child_Age);
    }
    
    return hr;
    
}


/**
 *
 */
void BurstEventProposal::prepareProposal( void )
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
void BurstEventProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void BurstEventProposal::undoProposal( void )
{
    
    // undo the proposal (only if succeeded)
    if ( failed == false)
    {
        BirthDeathBurstProcess *dist = static_cast<BirthDeathBurstProcess *>( &psi->getDistribution() );
        size_t node_index= stored_node->getIndex();
        
        if ( dist->isBurstSpeciation( node_index )  == true )
        {
            dist->setBurstSpeciation( node_index, false );
            stored_node->setAge( stored_age );
        }
        else
        {
            dist->setBurstSpeciation( node_index, true );
            stored_node->setAge( stored_age );
        }
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BurstEventProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == psi )
    {
        psi = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    
}


void BurstEventProposal::setProposalTuningParameter(double tp)
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
void BurstEventProposal::tune( double rate )
{
    
}

