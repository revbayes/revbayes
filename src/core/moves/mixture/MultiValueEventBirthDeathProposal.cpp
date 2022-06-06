#include "MultiValueEventBirthDeathProposal.h"

#include <cmath>

#include "AutocorrelatedEventDistribution.h"
#include "MultiValueEventDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MultiValueEvent.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MultiValueEventBirthDeathProposal::MultiValueEventBirthDeathProposal( StochasticNode<MultiValueEvent> *n ) : Proposal(),
    event_var( n )
{
    
    // tell the base class to add the node
    addNode( event_var );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MultiValueEventBirthDeathProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MultiValueEventBirthDeathProposal* MultiValueEventBirthDeathProposal::clone( void ) const
{
    
    return new MultiValueEventBirthDeathProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MultiValueEventBirthDeathProposal::getProposalName( void ) const
{
    static std::string name = "MultiValueEventBirthDeath";
    
    return name;
}


double MultiValueEventBirthDeathProposal::getProposalTuningParameter( void ) const
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
double MultiValueEventBirthDeathProposal::doProposal( void )
{
    
    double hr = 0.0;
    
    const MultiValueEventDistribution* dist_mve = dynamic_cast< const MultiValueEventDistribution*>( &event_var->getDistribution() );
    if ( dist_mve != NULL )
    {
        hr = doUncorrelatedProposal( dist_mve );
    }
    else
    {
        const AutocorrelatedEventDistribution* dist_ace = dynamic_cast< const AutocorrelatedEventDistribution*>( &event_var->getDistribution() );
        hr = doAutocorrelatedProposal( dist_ace );
    }
    
    
    return hr;
}


double MultiValueEventBirthDeathProposal::doAutocorrelatedProposal(const AutocorrelatedEventDistribution *d)
{
    const AutocorrelatedEventDistribution& dist_mve = *d;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MultiValueEvent &mve = event_var->getValue();
    long n_events = mve.getNumberOfEvents();
    
    double hr = 0.0;

    // We need to randomly pick a birth or death move
    // Otherwise we might give birth and die every time
    double u = rng->uniform01();
    
    if ( u > 0.5 || n_events == 0 )
    {
        // we pick a birth move
        was_birth = true;
        
        // increment the number of events
        mve.setNumberOfEvents( n_events + 1 );
        
        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        for (size_t i=0; i<priors.size(); ++i)
        {
            priors[i]->redrawValue();
            double new_val = priors[i]->getValue();
            mve.getValues(i).push_back( new_val );
            
            hr -= priors[i]->computeLnProbability();
        }
        
        if ( n_events == 0 )
        {
            hr += log(0.5);
        }
        
    }
    else
    {
        
        // we picked a death move
        was_birth = false;
        
        // decrement the number of events
        mve.setNumberOfEvents( n_events - 1 );
        
        // get the offsets
        const std::vector<long> &offset = dist_mve.getMinimumNumberOfEvents();
        
        // randomly pick an index
        size_t idx = floor( n_events * rng->uniform01() );
        stored_index = idx;
        
        // store the current values
        stored_values.clear();
        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        for (size_t i=0; i<mve.getNumberOfValues(); ++i)
        {
            size_t this_index = idx+offset[i];
            std::vector<double> &this_values = mve.getValues(i);
            double old_val = this_values[this_index];
            stored_values.push_back( old_val );
            
            priors[i]->setValue( new double(old_val) );
            hr += priors[i]->computeLnProbability();
            
            this_values.erase( this_values.begin()+this_index );
        }
        
        
        if ( n_events == 1 )
        {
            hr -= log(0.5);
        }
    }
    
    return hr;
}


double MultiValueEventBirthDeathProposal::doUncorrelatedProposal(const MultiValueEventDistribution *d)
{
    const MultiValueEventDistribution& dist_mve = *d;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    MultiValueEvent &mve = event_var->getValue();
    long n_events = mve.getNumberOfEvents();
    
    double hr = 0.0;

    // We need to randomly pick a birth or death move
    // Otherwise we might give birth and die every time
    double u = rng->uniform01();
    
    if ( u > 0.5 || n_events == 0 )
    {
        // we pick a birth move
        was_birth = true;
        
        // increment the number of events
        mve.setNumberOfEvents( n_events + 1 );
        
        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        for (size_t i=0; i<priors.size(); ++i)
        {
            priors[i]->redrawValue();
            double new_val = priors[i]->getValue();
            mve.getValues(i).push_back( new_val );
            
            hr -= priors[i]->computeLnProbability();
        }
        
        if ( n_events == 0 )
        {
            hr += log(0.5);
        }
        
    }
    else
    {
        
        // we picked a death move
        was_birth = false;
        
        // decrement the number of events
        mve.setNumberOfEvents( n_events - 1 );
        
        // get the offsets
        const std::vector<long> &offset = dist_mve.getMinimumNumberOfEvents();
        
        // randomly pick an index
        size_t idx = floor( n_events * rng->uniform01() );
        stored_index = idx;
        
        // store the current values
        stored_values.clear();
        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        for (size_t i=0; i<mve.getNumberOfValues(); ++i)
        {
            size_t this_index = idx+offset[i];
            std::vector<double> &this_values = mve.getValues(i);
            double old_val = this_values[this_index];
            stored_values.push_back( old_val );
            
            priors[i]->setValue( new double(old_val) );
            hr += priors[i]->computeLnProbability();
            
            this_values.erase( this_values.begin()+this_index );
        }
        
        
        if ( n_events == 1 )
        {
            hr -= log(0.5);
        }
    }
    
    return hr;
}


/**
 *
 */
void MultiValueEventBirthDeathProposal::prepareProposal( void )
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
void MultiValueEventBirthDeathProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MultiValueEventBirthDeathProposal::undoProposal( void )
{
    
    MultiValueEvent &mve = event_var->getValue();
    const MultiValueEventDistribution &dist_mve = static_cast< const MultiValueEventDistribution &>( event_var->getDistribution() );
    long n_events = mve.getNumberOfEvents();
    
    // undo the proposal
    if ( was_birth == true )
    {
        
        // decrement the number of events
        mve.setNumberOfEvents( n_events - 1 );
        
        // get the offsets
        const std::vector<long> &offset = dist_mve.getMinimumNumberOfEvents();
        
        // remove the proposed values
        for (size_t i=0; i<mve.getNumberOfValues(); ++i)
        {
            size_t this_index = n_events+offset[i]-1;
            std::vector<double> &this_values = mve.getValues(i);
            this_values.erase( this_values.begin()+this_index );
        }
        
    }
    else
    {
        // increment the number of events
        mve.setNumberOfEvents( n_events + 1 );
        
        // get the offsets
        const std::vector<long> &offset = dist_mve.getMinimumNumberOfEvents();

        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        for (size_t i=0; i<priors.size(); ++i)
        {
            size_t this_index = stored_index+offset[i];
            std::vector<double> &this_values = mve.getValues(i);
            this_values.insert ( this_values.begin()+this_index , stored_values[i] );
        }
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MultiValueEventBirthDeathProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == event_var )
    {
        event_var = static_cast<StochasticNode<MultiValueEvent>* >(newN) ;
    }
    
}


void MultiValueEventBirthDeathProposal::setProposalTuningParameter(double tp)
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
void MultiValueEventBirthDeathProposal::tune( double rate )
{
    
}

