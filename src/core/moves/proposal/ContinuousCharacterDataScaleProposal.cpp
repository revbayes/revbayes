#include "ContinuousCharacterDataScaleProposal.h"

#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ContinuousCharacterDataScaleProposal::ContinuousCharacterDataScaleProposal( StochasticNode<ContinuousCharacterData> *n, double l) : Proposal(),
variable( n ),
lambda( l ),
stored_value( 0.0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void ContinuousCharacterDataScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ContinuousCharacterDataScaleProposal* ContinuousCharacterDataScaleProposal::clone( void ) const
{
    
    return new ContinuousCharacterDataScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ContinuousCharacterDataScaleProposal::getProposalName( void ) const
{
    static std::string name = "MatrixRealSingleElementSlidingMove";
    
    return name;
}


double ContinuousCharacterDataScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * lambda = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double ContinuousCharacterDataScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    ContinuousCharacterData& v = variable->getValue();
    // choose an index
    index_a = size_t( rng->uniform01() * v.getNumberOfTaxa() );
    index_b = size_t( rng->uniform01() * v.getNumberOfCharacters() );
    
    // copy value
    stored_value = v[index_a][index_b];
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    v[index_a][index_b] *= scaling_factor;
    
    // compute the Hastings ratio
    double ln_hastings_ratio = log( scaling_factor );
    
    return ln_hastings_ratio;
}


/**
 *
 */
void ContinuousCharacterDataScaleProposal::prepareProposal( void )
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
void ContinuousCharacterDataScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void ContinuousCharacterDataScaleProposal::undoProposal( void )
{
    
    ContinuousCharacterData& v = variable->getValue();
    v[index_a][index_b] = stored_value;
    //    variable->clearTouchedElementIndices();
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void ContinuousCharacterDataScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<ContinuousCharacterData>* >(newN) ;
    
}


void ContinuousCharacterDataScaleProposal::setProposalTuningParameter(double tp)
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
void ContinuousCharacterDataScaleProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
    lambda = fmin(10000, lambda);
    
}

