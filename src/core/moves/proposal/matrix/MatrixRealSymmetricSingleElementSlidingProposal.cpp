#include "MatrixRealSymmetricSingleElementSlidingProposal.h"
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
MatrixRealSymmetricSingleElementSlidingProposal::MatrixRealSymmetricSingleElementSlidingProposal( StochasticNode<MatrixReal> *n, double l) : Proposal(),
variable( n ),
delta( l ),
storedValue( 0.0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MatrixRealSymmetricSingleElementSlidingProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MatrixRealSymmetricSingleElementSlidingProposal* MatrixRealSymmetricSingleElementSlidingProposal::clone( void ) const
{
    
    return new MatrixRealSymmetricSingleElementSlidingProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MatrixRealSymmetricSingleElementSlidingProposal::getProposalName( void ) const
{
    static std::string name = "MatrixRealSingleElementSlidingMove";
    
    return name;
}


double MatrixRealSymmetricSingleElementSlidingProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = delta * u
 * where delta is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MatrixRealSymmetricSingleElementSlidingProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MatrixReal& v = variable->getValue();
    // choose an index
    indexa = size_t( rng->uniform01() * v.getNumberOfRows() );
    indexb = size_t( rng->uniform01() * v.getNumberOfColumns() );
    
    // copy value
    storedValue = v[indexa][indexb];
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double slidingFactor = delta * ( u - 0.5 );
    v[indexa][indexb] += slidingFactor;
    
    if (indexa != indexb)
    {
        v[indexb][indexa] += slidingFactor;
    }

    variable->addTouchedElementIndex(indexa*v.getNumberOfRows() + indexb);
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
}


/**
 *
 */
void MatrixRealSymmetricSingleElementSlidingProposal::prepareProposal( void )
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
void MatrixRealSymmetricSingleElementSlidingProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void MatrixRealSymmetricSingleElementSlidingProposal::undoProposal( void )
{
    
    MatrixReal& v = variable->getValue();
    v[indexa][indexb] = storedValue;
    v[indexb][indexa] = storedValue;
    variable->clearTouchedElementIndices();
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MatrixRealSymmetricSingleElementSlidingProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    
}


void MatrixRealSymmetricSingleElementSlidingProposal::setProposalTuningParameter(double tp)
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
void MatrixRealSymmetricSingleElementSlidingProposal::tune( double rate )
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

