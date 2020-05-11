#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionNormal.h"
#include "ContinuousEventScaleProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "CharacterEventContinuous.h"
#include "CharacterHistoryContinuous.h"
#include "ConditionedBirthDeathShiftProcessContinuous.h"
#include "Move.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ContinuousEventScaleProposal::ContinuousEventScaleProposal( StochasticNode<Tree> *n, double l ) : Proposal(),
    variable( n ),
    lambda( l ),
    stored_value()
{
    // tell the base class to add the node
    addNode( variable );
    
    distribution = dynamic_cast< ConditionedBirthDeathShiftProcessContinuous* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for continuous-event shift move.");
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ContinuousEventScaleProposal::cleanProposal( void )
{
    
    // reset the stored value
    stored_value = NULL;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ContinuousEventScaleProposal* ContinuousEventScaleProposal::clone( void ) const
{
    
    return new ContinuousEventScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ContinuousEventScaleProposal::getProposalName( void ) const
{
    static std::string name = "DiscreteEventCategoryRandomWalk";
    
    return name;
}


double ContinuousEventScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double ContinuousEventScaleProposal::doProposal( void )
{
    
    // define the Hastings ratio
    double ln_hastings_ratio = 0.0;
    
    CharacterHistoryContinuous &history = distribution->getCharacterHistory();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    std::vector<size_t> rate_indices;
    if ( distribution->isSpeciationRateConstant() == false )
    {
        rate_indices.push_back( 0 );
    }
    if ( distribution->isExtinctionRateConstant() == false )
    {
        rate_indices.push_back( 1 );
    }
    
    size_t num_events = history.getNumberEvents();
    failed = (num_events == 0 && rate_indices.size() > 0);
    
    if ( failed == false )
    {
        // pick a random event
        size_t branch_index = 0;
        CharacterEventContinuous *event = history.pickRandomEvent( branch_index );

        // pick a index (either speciation or extinction)
        size_t tmp = size_t( floor( rng->uniform01() * rate_indices.size() ) );
        size_t rate_index = rate_indices[tmp];

        // store the index
        stored_index = rate_index;
        // store the event
        stored_value = event;
        // get the current index
        stored_rate = event->getState(stored_index);
    
        // draw a new value
        double u = rng->uniform01();
        double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
        double new_value = stored_rate * scaling_factor;
    
        // compute the Hastings ratio
        ln_hastings_ratio += log( scaling_factor );

        // finally, set the new value
        event->setState( new_value, stored_index );
        
    }
    else
    {
        move->decrementTriedCounter();
    }
    
    return ln_hastings_ratio;
}


/**
 *
 */
void ContinuousEventScaleProposal::prepareProposal( void )
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
void ContinuousEventScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void ContinuousEventScaleProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        stored_value->setState( stored_rate, stored_index );
        stored_value = NULL;
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ContinuousEventScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< ConditionedBirthDeathShiftProcessContinuous* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for continuous-event scale move.");
    }
}


void ContinuousEventScaleProposal::setProposalTuningParameter(double tp)
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
void ContinuousEventScaleProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
}


