#include "ContinuousCharacterDataSlideProposal.h"
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
ContinuousCharacterDataSlideProposal::ContinuousCharacterDataSlideProposal( StochasticNode<ContinuousCharacterData> *n, double l) : Proposal(),
    variable( n ),
    delta( l ),
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
void ContinuousCharacterDataSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ContinuousCharacterDataSlideProposal* ContinuousCharacterDataSlideProposal::clone( void ) const
{
    
    return new ContinuousCharacterDataSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& ContinuousCharacterDataSlideProposal::getProposalName( void ) const
{
    static std::string name = "ContinuousCharacterDataSlideMove";
    
    return name;
}


double ContinuousCharacterDataSlideProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double ContinuousCharacterDataSlideProposal::doProposal( void )
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
    double sliding_factor = delta * ( u - 0.5 );
    v[index_a][index_b] += sliding_factor;
    
    
//    variable->addTouchedElementIndex(indexa);
//    variable->addTouchedElementIndex(indexb);
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
}


/**
 *
 */
void ContinuousCharacterDataSlideProposal::prepareProposal( void )
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
void ContinuousCharacterDataSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void ContinuousCharacterDataSlideProposal::undoProposal( void )
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
void ContinuousCharacterDataSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<ContinuousCharacterData>* >(newN) ;
    
}


void ContinuousCharacterDataSlideProposal::setProposalTuningParameter(double tp)
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
void ContinuousCharacterDataSlideProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta /= (2.0 - rate/0.44 );
    }
    
    delta = fmin(10000, delta);
    
}

