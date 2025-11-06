#include "MultiValueEventSlideProposal.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "Cloneable.h"
#include "Move.h"
#include "MultiValueEvent.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MultiValueEventSlideProposal::MultiValueEventSlideProposal( StochasticNode<MultiValueEvent> *n, const std::string &vn, double l ) : Proposal(),
    event_var( n ),
    value_name( vn ),
    failed( false ),
    lambda( l )
{
    
    // tell the base class to add the node
    addNode( event_var );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MultiValueEventSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MultiValueEventSlideProposal* MultiValueEventSlideProposal::clone( void ) const
{
    
    return new MultiValueEventSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MultiValueEventSlideProposal::getProposalName( void ) const
{
    static std::string name = "MultiValueEventSlide";
    
    return name;
}


double MultiValueEventSlideProposal::getProposalTuningParameter( void ) const
{
    
    return lambda;
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
double MultiValueEventSlideProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MultiValueEvent &mve = event_var->getValue();
    std::vector<double> &values = mve.getValues( value_name );
    
    double hr = RbConstants::Double::neginf;
    if ( values.size() > 0 )
    {
        failed = false;
        
        // choose an index
        size_t index = size_t(rng->uniform01() * values.size());
        
        double &val = values[index];
        
        // copy value
        stored_value = val;
        stored_index = index;
        
        // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
        double u       = rng->uniform01();
        double delta   = ( lambda * ( u - 0.5 ) );
        double new_val = val + delta;
        
        val = new_val;
        
        hr = 0.0;
    }
    else
    {
        failed = true;
        // we also need to decrement the counters for our tuning
        move->decrementTriedCounter();
    }
    
    return hr;
}


/**
 *
 */
void MultiValueEventSlideProposal::prepareProposal( void )
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
void MultiValueEventSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void MultiValueEventSlideProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        MultiValueEvent &mve = event_var->getValue();
        std::vector<double> &values = mve.getValues( value_name );
        values[stored_index] = stored_value;
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MultiValueEventSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == event_var )
    {
        event_var = static_cast<StochasticNode<MultiValueEvent>* >(newN) ;
    }
    
}


void MultiValueEventSlideProposal::setProposalTuningParameter(double tp)
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
void MultiValueEventSlideProposal::tune( double rate )
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
    
    lambda = fmin(10000, lambda);
    
}
