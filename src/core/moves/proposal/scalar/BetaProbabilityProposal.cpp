#include "BetaProbabilityProposal.h"
#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbStatisticsHelper.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>


using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BetaProbabilityProposal::BetaProbabilityProposal( StochasticNode<double> *n, double d, double o ) : SimpleProposal<double>( n ),
    stored_value( 0.0 ),
    delta( d ),
    offset( o )
{
    // tell the base class to add the node
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BetaProbabilityProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BetaProbabilityProposal* BetaProbabilityProposal::clone( void ) const
{
    
    return new BetaProbabilityProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BetaProbabilityProposal::getProposalName( void ) const
{
    static std::string name = "BetaProbability";
    
    return name;
}


double BetaProbabilityProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * A Beta-Probability proposal randomly changes some values of a Probability, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new Probability
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The Probability is then renormalized.
 *
 * \return The hastings ratio.
 */
double BetaProbabilityProposal::propose( double &value )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // store the value
    stored_value = value;
    
    // get the current value
    double current_value = value;
    
    // draw new rates and compute the hastings ratio at the same time
    double a = delta * current_value + offset;
    double b = (1.0-current_value) * delta + offset;
    double new_value = RbStatistics::Beta::rv(a, b, *rng);
    
    // set the value
    value = new_value;
    
    double ln_Hastins_ratio = 0;
    try
    {
        // compute the Hastings ratio
        double forward = RbStatistics::Beta::lnPdf(a, b, new_value);
        double new_a = delta * new_value + offset;
        double new_b = (1.0-new_value) * delta + offset;
        double backward = RbStatistics::Beta::lnPdf(new_a, new_b, current_value);
        
        ln_Hastins_ratio = backward - forward;
    }
	catch (RbException &e)
    {
        ln_Hastins_ratio = RbConstants::Double::neginf;
    }
    
    return ln_Hastins_ratio;
}


/**
 *
 */
void BetaProbabilityProposal::prepareProposal( void )
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
void BetaProbabilityProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void BetaProbabilityProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( new double(stored_value) );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BetaProbabilityProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<double>* >(newN) ;
    
}


void BetaProbabilityProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void BetaProbabilityProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta /= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta *= (2.0 - rate/0.44 );
    }
    
}

