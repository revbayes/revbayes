#include <cmath>
#include <iostream>

#include "DistributionBeta.h"
#include "RandomDiveProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "ContinuousStochasticNode.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RandomDiveProposal::RandomDiveProposal( StochasticNode<double> *n, double l, double p) : Proposal(p),
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
void RandomDiveProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RandomDiveProposal* RandomDiveProposal::clone( void ) const
{
    
    return new RandomDiveProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RandomDiveProposal::getProposalName( void ) const
{
    static std::string name = "RandomDive";
    
    return name;
}


double RandomDiveProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A Random Dive proposal draws a multiplier, epsilon, from a distribution on (-1,1) - {0} distribution
 * and multiplies the current value by either epsilon or 1/epsilon with equal probability,
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 * Here we use a beta distribution (stretched/shifted to fit on (-1,1)) so that we can tune the proposal.
 * For more details, see Dutta, Somak. "Multiplicative random walk Metropolis-Hastings on the real line."
 *
 * \return The hastings ratio.
 */
double RandomDiveProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = variable->getValue();
    
    // copy value
    storedValue = val;
    
//    // Draw from a stretched and translated Beta(lambda,lambda) - {0} distribution that does not include -1 or 1
//    double epsilon  = ( 2 * RbStatistics::Beta::rv(lambda, lambda, *rng) ) - 1.0 ;
    
    // Draw from a mirrored Beta(lambda,1) - {0} distribution that does not include -1 or 1
    
    // Draw
    double epsilon  = RbStatistics::Beta::rv(lambda, 1.0, *rng);

    while ( epsilon == 0.0 || epsilon == 1.0 )
    {
        epsilon  = RbStatistics::Beta::rv(lambda, 1.0, *rng);
    }
    
    // Mirror
    if ( rng->uniform01() < 0.5 )
    {
        epsilon *= -1.0;
    }
    
    double u = rng->uniform01();
    
    if (u >= 0.5)
    {
        epsilon = 1/epsilon;
    }
    
    double newVal = val * epsilon;
    
    val = newVal;
    
    return log(fabs(epsilon));
}


/**
 *
 */
void RandomDiveProposal::prepareProposal( void )
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
void RandomDiveProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
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
void RandomDiveProposal::undoProposal( void )
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
void RandomDiveProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< ContinuousStochasticNode* >(newN) ;
    
}


void RandomDiveProposal::setProposalTuningParameter(double tp)
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
void RandomDiveProposal::tune( double rate )
{
    // As the proposal distribution gets more concentrated around 0, the proposals get bolder
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= 0.95;
    }
    else
    {
        lambda *= 1.05;
    }
    
    lambda = fmin(10000, lambda);
    lambda = fmax(0.0001, lambda);
    
}

