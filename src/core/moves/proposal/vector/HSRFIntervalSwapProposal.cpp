#include "HSRFIntervalSwapProposal.h"

#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
HSRFIntervalSwapProposal::HSRFIntervalSwapProposal( std::vector< StochasticNode<double> *> d, std::vector< StochasticNode<double> *> s) : Proposal(),
    delta( d ),
    sigma( s ),
    stored_index(  )
{
    // tell the base class to add the node
    for (std::vector< StochasticNode<double> *>::const_iterator it = delta.begin(); it != delta.end(); it++)
    {
        addNode( *it );
    }
    for (std::vector< StochasticNode<double> *>::const_iterator it = sigma.begin(); it != sigma.end(); it++)
    {
        addNode( *it );
    }

    if ( delta.size() != sigma.size() )
    {
      throw(RbException("HSRFIntervalSwapProposal requires that both vectors are of the same length."));
    }

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void HSRFIntervalSwapProposal::cleanProposal( void )
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
HSRFIntervalSwapProposal* HSRFIntervalSwapProposal::clone( void ) const
{

    return new HSRFIntervalSwapProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& HSRFIntervalSwapProposal::getProposalName( void ) const
{
    static std::string name = "HSRFIntervalSwap";

    return name;
}


double HSRFIntervalSwapProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
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
double HSRFIntervalSwapProposal::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    // Do we first swap 1,2 or 2,3? 50-50 for either, so move is symmetric
    double u = rng->uniform01() * 2;
    size_t start = size_t(floor(u));
    size_t end_plus_1 = size_t(start + 2 * (floor((delta.size() - start)/2)-1) + 1);

    // store the starting index so we can revert the move
    stored_index = start;
    size_t i = start;
    for (size_t i=start; i<end_plus_1; ++i) {

      // copy value
      double d_i = delta[i]->getValue();
      double s_i = sigma[i]->getValue();

      // Swap
      delta[i]->getValue() = delta[i+1]->getValue();
      sigma[i]->getValue() = sigma[i+1]->getValue();

      delta[i+1]->getValue() = d_i;
      sigma[i+1]->getValue() = s_i;

      ++i;
    }


    return 0.0;

}


/**
 * Prepare the proposal, e.g., pick the element that we want to change.
 * Here we do not need to do any preparation.
 */
void HSRFIntervalSwapProposal::prepareProposal( void )
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
void HSRFIntervalSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void HSRFIntervalSwapProposal::undoProposal( void )
{

    size_t start = stored_index;
    size_t end_plus_1 = size_t(start + 2 * (floor((delta.size() - start)/2)-1) + 1);

    // Undoing the move is just running the move exactly the same way (1,2),(3,4)->(2,1),(4,3)->(1,2),(3,4)
    for (size_t i=start; i<end_plus_1; ++i) {

      // copy value
      double d_i = delta[i]->getValue();
      double s_i = sigma[i]->getValue();

      // Swap
      delta[i]->getValue() = delta[i+1]->getValue();
      sigma[i]->getValue() = sigma[i+1]->getValue();

      delta[i+1]->getValue() = d_i;
      sigma[i+1]->getValue() = s_i;

      ++i;
    }

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void HSRFIntervalSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    for (size_t i = 0; i < delta.size(); ++i)
    {
        if ( delta[i] == oldN )
        {
            delta[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

    for (size_t i = 0; i < delta.size(); ++i)
    {
        if ( sigma[i] == oldN )
        {
            sigma[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

}


void HSRFIntervalSwapProposal::setProposalTuningParameter(double tp)
{

}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void HSRFIntervalSwapProposal::tune( double rate )
{

}
