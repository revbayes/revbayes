#include "DistributionBeta.h"
#include "HalfRandomDiveProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
HalfRandomDiveProposal::HalfRandomDiveProposal( StochasticNode<double> *n, double l, double p) : Proposal(p),
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
void HalfRandomDiveProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
HalfRandomDiveProposal* HalfRandomDiveProposal::clone( void ) const
{
    
    return new HalfRandomDiveProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& HalfRandomDiveProposal::getProposalName( void ) const
{
    static std::string name = "HalfRandomDive";
    
    return name;
}


double HalfRandomDiveProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A Random Dive proposal on the that does not cross 0.
 *
 * \return The hastings ratio.
 */
double HalfRandomDiveProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = variable->getValue();
    
    // copy value
    storedValue = val;
    
    // Draw from a Beta(lambda,1) - {0} distribution that does not 1
    double epsilon  = RbStatistics::Beta::rv(lambda, 1.0, *rng);
    while ( epsilon == 0.0 || epsilon == 1.0 )
    {
        epsilon  = RbStatistics::Beta::rv(lambda, 1.0, *rng);
    }
    
    double u = rng->uniform01();
    
    if (u >= 0.5)
    {
        epsilon = 1/epsilon;
    }
    
    double hr = log(fabs(epsilon));
    
    val = val * epsilon;
    
    return hr;
}


/**
 *
 */
void HalfRandomDiveProposal::prepareProposal( void )
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
void HalfRandomDiveProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void HalfRandomDiveProposal::undoProposal( void )
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
void HalfRandomDiveProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< ContinuousStochasticNode* >(newN) ;
    
}


void HalfRandomDiveProposal::setProposalTuningParameter(double tp)
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
void HalfRandomDiveProposal::tune( double rate )
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

