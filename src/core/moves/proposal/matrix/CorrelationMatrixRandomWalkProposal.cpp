#include "CorrelationMatrixRandomWalkProposal.h"

#include <cstddef>
#include <iostream>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
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
CorrelationMatrixRandomWalkProposal::CorrelationMatrixRandomWalkProposal( StochasticNode<MatrixReal> *n, double s, double p) : Proposal(p),
    variable( n ),
    sigma( s ),
    stored_matrix()
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CorrelationMatrixRandomWalkProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixRandomWalkProposal* CorrelationMatrixRandomWalkProposal::clone( void ) const
{
    
    return new CorrelationMatrixRandomWalkProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixRandomWalkProposal::getProposalName( void ) const
{
    static std::string name = "CorrelationMatrixRandomWalk";
    
    return name;
}


double CorrelationMatrixRandomWalkProposal::getProposalTuningParameter( void ) const
{
    return sigma;
}


/**
 * Perform the proposal.
 *
 * A Beta proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementBetas the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double CorrelationMatrixRandomWalkProposal::doProposal( void )
{
    
    stored_matrix = variable->getValue();
    size_t dim = stored_matrix.size();
    
    MatrixReal& current_matrix = variable->getValue();
    for (size_t i = 0; i < (dim - 1); ++i)
    {
        for (size_t j = i+1; j < dim; ++j)
        {
            current_matrix[i][j] = RbStatistics::Normal::rv(stored_matrix[i][j], sigma, *GLOBAL_RNG);
            current_matrix[j][i] = current_matrix[i][j];
            variable->addTouchedElementIndex(i);
            variable->addTouchedElementIndex(j);
        }
    }
    variable->touch(true);
    
    return 0.0;
    
}


/**
 *
 */
void CorrelationMatrixRandomWalkProposal::prepareProposal( void )
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
void CorrelationMatrixRandomWalkProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "sigma = ";
    if (name_only == false)
    {
        o << sigma;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void CorrelationMatrixRandomWalkProposal::undoProposal( void )
{
    
    MatrixReal& current_matrix = variable->getValue();
    size_t size = current_matrix.getNumberOfRows();
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = i; j < size; ++j){
            current_matrix[i][j] = stored_matrix[i][j];
            current_matrix[j][i] = stored_matrix[i][j];
        }
    }
    variable->clearTouchedElementIndices();
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void CorrelationMatrixRandomWalkProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    
}


void CorrelationMatrixRandomWalkProposal::setProposalTuningParameter(double tp)
{
    sigma = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void CorrelationMatrixRandomWalkProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    
    if ( rate > p )
    {
        sigma *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        sigma /= (2.0 - rate/p);
    }
    
}

