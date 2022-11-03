#include "CorrelationMatrixElementSwapProposal.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
CorrelationMatrixElementSwapProposal::CorrelationMatrixElementSwapProposal( StochasticNode<MatrixReal> *n) : Proposal(),
    variable( n ),
    storedValue_1( 0.0 ),
    storedValue_2( 0.0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CorrelationMatrixElementSwapProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixElementSwapProposal* CorrelationMatrixElementSwapProposal::clone( void ) const
{
    
    return new CorrelationMatrixElementSwapProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixElementSwapProposal::getProposalName( void ) const
{
    static std::string name = "CorrelationMatrixElementSwapMove";
    
    return name;
}


/**
 * Perform the proposal.
 *
 * A Swap proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSwaps the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double CorrelationMatrixElementSwapProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MatrixReal& v = variable->getValue();
    
    // choose an index for the first element
    indexa_1 = size_t( rng->uniform01() * v.getNumberOfRows() );
    indexb_1 = size_t( rng->uniform01() * v.getNumberOfColumns() );
    
    // make sure we don't get a diagonal
    while (indexb_1 == indexa_1)
    {
        indexb_1 = size_t( rng->uniform01() * v.getNumberOfColumns() );
    }

    // choose an index for the second element
    indexa_2 = size_t( rng->uniform01() * v.getNumberOfRows() );
    indexb_2 = size_t( rng->uniform01() * v.getNumberOfColumns() );
    
    // make sure we don't get a diagonal
    while (indexb_2 == indexa_2)
    {
        indexb_2 = size_t( rng->uniform01() * v.getNumberOfColumns() );
    }
    
    // at this point, we may choose the same element across the diagonal...
    
    // copy the current values
    storedValue_1 = v[indexa_1][indexb_1];
    storedValue_2 = v[indexa_2][indexb_2];

    // swap the values
    double current_value_1 = storedValue_1;
    double current_value_2 = storedValue_2;
    
    v[indexa_1][indexb_1] = current_value_2;
    v[indexb_1][indexa_1] = current_value_2;
    
    v[indexa_2][indexb_2] = current_value_1;
    v[indexb_2][indexa_2] = current_value_1;
    
    // touch the elements of the matrix
    variable->addTouchedElementIndex(indexa_1);
    variable->addTouchedElementIndex(indexb_1);
    variable->addTouchedElementIndex(indexa_2);
    variable->addTouchedElementIndex(indexb_2);

    // this proposal is symmetrical
    return 0.0;
    
}


/**
 *
 */
void CorrelationMatrixElementSwapProposal::prepareProposal( void )
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
void CorrelationMatrixElementSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void CorrelationMatrixElementSwapProposal::undoProposal( void )
{
    
    MatrixReal& v = variable->getValue();
    
    
    v[indexa_1][indexb_1] = storedValue_1;
    v[indexb_1][indexa_1] = storedValue_1;

    v[indexa_2][indexb_2] = storedValue_2;
    v[indexb_2][indexa_2] = storedValue_2;

    variable->clearTouchedElementIndices();
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void CorrelationMatrixElementSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    
}
