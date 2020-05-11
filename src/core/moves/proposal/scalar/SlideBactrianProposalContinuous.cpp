#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "SlideBactrianProposalContinuous.h"
#include "ContinuousStochasticNode.h"
#include "Proposal.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SlideBactrianProposalContinuous::SlideBactrianProposalContinuous( ContinuousStochasticNode *n, double l, double r) : Proposal(r),
    variable( n ),
    storedValue( 0.0 ),
    lambda( l )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SlideBactrianProposalContinuous::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SlideBactrianProposalContinuous* SlideBactrianProposalContinuous::clone( void ) const
{
    
    return new SlideBactrianProposalContinuous( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SlideBactrianProposalContinuous::getProposalName( void ) const
{
    static std::string name = "SlidingBactrian";
    
    return name;
}


double SlideBactrianProposalContinuous::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and slides the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double SlideBactrianProposalContinuous::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = variable->getValue();
    
    // copy value
    storedValue = val;
    
    double min = variable->getMin();
    double max = variable->getMax();
        
    //double size = max - min;
    
    double u = rng->uniform01();
    double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * RbConstants::BACT_SD;
    
    // See Yang and Rodriguez (2013) SI eqns 19 and 20
    // Currently hard-coding m = 0.95
    if ( u < 0.5 ) {
        delta += 0.95;
    } else {
        delta -= 0.95;
    }

    double newVal = val + (lambda * delta);
    
    /* reflect the new value */
    do {
        if ( newVal < min )
        {
            newVal = 2.0 * min - newVal;
        }
        else if ( newVal > max )
        {
            newVal = 2.0 * max - newVal;
        }
    } while ( newVal < min || newVal > max );
    
    // FIXME: not the most efficient way of handling multiple reflections :-P
    
    val = newVal;
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
}


/**
 *
 */
void SlideBactrianProposalContinuous::prepareProposal( void )
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
void SlideBactrianProposalContinuous::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "sigma = ";
    if (name_only == false)
    {
        o << lambda;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SlideBactrianProposalContinuous::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( new double(storedValue) );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SlideBactrianProposalContinuous::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< ContinuousStochasticNode* >(newN) ;
    
}


void SlideBactrianProposalContinuous::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SlideBactrianProposalContinuous::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        lambda /= (2.0 - rate/p);
    }
    
    if ( lambda > 1 ) {
        lambda = fmin(10000, lambda);
    } else {
        lambda = fmax(1/10000, lambda);
    }
    
}

