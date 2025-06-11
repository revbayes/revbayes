#include "MultiValueEventBirthDeathProposal.h"

#include <cmath>

#include "AutocorrelatedEventDistribution.h"
#include "DistributionNormal.h"
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
MultiValueEventBirthDeathProposal::MultiValueEventBirthDeathProposal( StochasticNode<MultiValueEvent> *n, bool use_ac ) : Proposal(),
    event_var( n ),
    use_autocorrelated_proposal( use_ac )
{
    
    // tell the base class to add the node
    addNode( event_var );
    
    size_t num_vals = n->getValue().getNumberOfValues();
    ac_proposal_sd = std::vector<double>(num_vals, 0.1);
    
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
    std::int64_t n_events = mve.getNumberOfEvents();
    
    double hr = 0.0;

    // We need to randomly pick a birth or death move
    // Otherwise we might give birth and die every time
    double u = rng->uniform01();
    
    stored_sorting_index = -1;
    
    if ( u > 0.5 || n_events == 0 )
    {
        // we pick a birth move
        was_birth = true;
        
        // increment the number of events
        mve.setNumberOfEvents( n_events + 1 );
        
        // get the offsets
        const std::vector<std::int64_t> &offset = dist_mve.getMinimumNumberOfEvents();
        
        std::vector< TypedDistribution<double> * > priors = dist_mve.getValuePriors();
        
        for (size_t i=0; i<priors.size(); ++i)
        {
            std::vector<double>& vals = mve.getValues(i);
            
            if ( dist_mve.isSorted( i ) == true )
            {
                size_t num_vals = vals.size();
                int left  = -1;
                int right = int(num_vals);
                                
                bool success = false;
                double new_val = 0.0;
                while ( success == false )
                {
                    if ( num_vals > 0 )
                    {
                        left = int( (num_vals+1) * rng->uniform01()) - 1;
                        right = left + 1;
                    }
                    
                    success = true;
                    priors[i]->redrawValue();
                    new_val = priors[i]->getValue();
                    if ( left > -1 )
                    {
                        success &= ( new_val > vals[left] );
                    }
                    if ( right < num_vals )
                    {
                        success &= ( new_val < vals[right] );
                    }
                }
                
                stored_sorting_index = left+1;
                vals.insert(vals.begin()+stored_sorting_index, new_val);
                
                hr -= priors[i]->computeLnProbability();
                
            }
            else
            {

                
                size_t insert_index = offset[i];
                if ( stored_sorting_index != -1 )
                {
                    insert_index += stored_sorting_index;
                }

                double new_val = -1.0;
                if ( dist_mve.isAutocorrelated( i ) == false || use_autocorrelated_proposal == false )
                {
                    priors[i]->redrawValue();
                    new_val = priors[i]->getValue();
                    hr -= priors[i]->computeLnProbability();
                }
                else
                {
                    double predecessor_value = vals[insert_index-1];
                    new_val = RbStatistics::Normal::rv(predecessor_value, ac_proposal_sd[i], *GLOBAL_RNG);
                    hr -= RbStatistics::Normal::lnPdf(predecessor_value, ac_proposal_sd[i], new_val);
                }

                vals.insert(vals.begin()+insert_index, new_val);
            }
            
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
        const std::vector<std::int64_t> &offset = dist_mve.getMinimumNumberOfEvents();
        
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

            if ( dist_mve.isAutocorrelated( i ) == false || use_autocorrelated_proposal == false )
            {
                priors[i]->setValue( new double(old_val) );
                hr += priors[i]->computeLnProbability();
            }
            else
            {
                double predecessor_value = this_values[this_index-1];
                hr += RbStatistics::Normal::lnPdf(predecessor_value, ac_proposal_sd[i], old_val);
            }
            
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
    std::int64_t n_events = mve.getNumberOfEvents();
    
    double hr = 0.0;

    // We need to randomly pick a birth or death move
    // Otherwise we might give birth and die every time
    double u = rng->uniform01();
    
    stored_sorting_index = -1;
    
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
        const std::vector<std::int64_t> &offset = dist_mve.getMinimumNumberOfEvents();
        
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
    
    o << "sd = ";
    if (name_only == false)
    {
        o << ac_proposal_sd[0];
    }
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
    std::int64_t n_events = mve.getNumberOfEvents();
    
    // undo the proposal
    if ( was_birth == true )
    {
        
        // decrement the number of events
        mve.setNumberOfEvents( n_events - 1 );
        
        // get the offsets
        const std::vector<std::int64_t> &offset = dist_mve.getMinimumNumberOfEvents();
        
        // remove the proposed values
        for (size_t i=0; i<mve.getNumberOfValues(); ++i)
        {
            size_t this_index = n_events+offset[i]-1;
            std::vector<double> &this_values = mve.getValues(i);
            if ( stored_sorting_index != -1 )
            {
                this_values.erase( this_values.begin()+stored_sorting_index+offset[i] );
            }
            else
            {
                this_values.erase( this_values.begin()+this_index );
            }
        }
        
    }
    else
    {
        // increment the number of events
        mve.setNumberOfEvents( n_events + 1 );
        
        // get the offsets
        const std::vector<std::int64_t> &offset = dist_mve.getMinimumNumberOfEvents();

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
    
    // Sebastian: auto-tuning doesn't seem to work for this move.
//    double p = this->targetAcceptanceRate;
//    if ( rate > p )
//    {
//        for ( size_t i=0; i<ac_proposal_sd.size(); ++i )
//        {
//            ac_proposal_sd[i] *= (1.0 + ((rate-p)/(1.0 - p)) );
//            ac_proposal_sd[i] = ( ac_proposal_sd[i] > 100.0 ? 100.0 : ac_proposal_sd[i] );
//        }
//
//    }
//    else
//    {
//        for ( size_t i=0; i<ac_proposal_sd.size(); ++i )
//        {
//            ac_proposal_sd[i] /= (2.0 - rate/p);
//            ac_proposal_sd[i] = ( ac_proposal_sd[i] < 0.001 ? 0.001 : ac_proposal_sd[i] );
//        }
//    }
}
