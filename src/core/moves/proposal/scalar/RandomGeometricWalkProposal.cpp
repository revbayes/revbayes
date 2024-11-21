#include <iostream>
#include <cstdint>

#include "DistributionGeometric.h"
#include "RandomGeometricWalkProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RandomGeometricWalkProposal::RandomGeometricWalkProposal( StochasticNode<std::int64_t> *n, double a ) : Proposal(),
    variable( n ),
    stored_value( 0 ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RandomGeometricWalkProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RandomGeometricWalkProposal* RandomGeometricWalkProposal::clone( void ) const
{
    
    return new RandomGeometricWalkProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RandomGeometricWalkProposal::getProposalName( void ) const
{
    static std::string name = "RandomGeometricWalk";
    
    return name;
}


double RandomGeometricWalkProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (0.0,1.0)
 * and moves one step up if u > 0.5 and one step down otherwise.
 *
 * \return The hastings ratio.
 */
double RandomGeometricWalkProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    std::int64_t &val = variable->getValue();
    
    // copy value
    stored_value = val;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    if ( u >= 0.5 )
    {
        val += RbStatistics::Geometric::rv(alpha, *rng);
    }
    else
    {
        val -= RbStatistics::Geometric::rv(alpha, *rng);
    }
    
    return 0.0;
}


/**
 *
 */
void RandomGeometricWalkProposal::prepareProposal( void )
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
void RandomGeometricWalkProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RandomGeometricWalkProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->getValue() = stored_value;
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RandomGeometricWalkProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<std::int64_t>* >(newN) ;
    
}


void RandomGeometricWalkProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void RandomGeometricWalkProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        alpha -= ( alpha/2.0 * ((rate-0.44)/0.56) );
    }
    else
    {
        alpha += ( (1.0-alpha)/2.0 * rate/0.45 );
    }
    
}

