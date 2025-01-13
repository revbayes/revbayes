#include "ScaleProposal.h"

#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ScaleProposal::ScaleProposal( StochasticNode<double> *n, double l, double p) : SimpleProposal<double>( n, p ),
    storedValue( 0.0 ),
    lambda( l )
{
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal. 
 */
ScaleProposal* ScaleProposal::clone( void ) const 
{
    
    return new ScaleProposal( *this );
}


/**
 * Get Proposals' name of object 
 *
 * \return The Proposals' name.
 */
const std::string& ScaleProposal::getProposalName( void ) const 
{
    static std::string name = "Scaling";
    
    return name;
}


double ScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/** 
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double ScaleProposal::propose( double &val )
{
    
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
        
    // copy value
    storedValue = val;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    val *= scaling_factor;
    
    // compute the Hastings ratio
    double ln_hastings_ratio = log( scaling_factor );
    
    return ln_hastings_ratio;
}


/**
 *
 */
void ScaleProposal::prepareProposal( void )
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
void ScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "lambda = ";
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
void ScaleProposal::undoProposal( void ) 
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
void ScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<double>* >(newN) ;
    
}


void ScaleProposal::setProposalTuningParameter(double tp)
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
void ScaleProposal::tune( double rate ) 
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
    
    if(lambda > 1) lambda = fmin(1000, lambda);
    else lambda = fmax(0.001, lambda);
}

