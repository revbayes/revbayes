#include "CorrelationMatrixPartialElementBetaProposal.h"

#include <cstddef>
#include <cmath>
#include <iostream>

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbConstants.h"
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
CorrelationMatrixPartialElementBetaProposal::CorrelationMatrixPartialElementBetaProposal( StochasticNode<MatrixReal> *n, double a) : Proposal(),
    variable( n ),
    alpha( a ),
    stored_matrix( n->getValue().getDim() )
{
    // tell the base class to add the node
    addNode( variable );
    
    // initialize the stored value
//    stored_matrix = variable->getValue();
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CorrelationMatrixPartialElementBetaProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixPartialElementBetaProposal* CorrelationMatrixPartialElementBetaProposal::clone( void ) const
{
    
    return new CorrelationMatrixPartialElementBetaProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixPartialElementBetaProposal::getProposalName( void ) const
{
    static std::string name = "MatrixRealSingleElementBetaMove";
    
    return name;
}


double CorrelationMatrixPartialElementBetaProposal::getProposalTuningParameter( void ) const
{
    return alpha;
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
double CorrelationMatrixPartialElementBetaProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // first, store the old matrix
    stored_matrix = variable->getValue();
    
    size_t dim = stored_matrix.getDim();
    
    // now, get the inverse of the matrix
    variable->getValue().setCholesky(true);
    MatrixReal inverse_matrix = variable->getValue().computeInverse();
    
    // compute the correlation matrix
    MatrixReal partial_correlations(dim, dim);

    for (size_t i = 0; i < dim; ++i)
    {
        partial_correlations[i][i] = 1.0;
        for (size_t j = i + 1; j < dim; ++j)
        {
            partial_correlations[i][j] = -inverse_matrix[i][j] / std::pow( inverse_matrix[i][i] * inverse_matrix[j][j], 0.5);
            partial_correlations[j][i] = -inverse_matrix[j][i] / std::pow( inverse_matrix[i][i] * inverse_matrix[j][j], 0.5);
        }
    }

    // compute the Jacobian of the transformation from a correlation matrix
    // to a partial correlation matrix
    double ln_jacobian = 0.0;
    for (size_t k = 0; k < dim - 1; ++k)
    {
        for (size_t i = k + 1; i < dim; ++i)
        {
            ln_jacobian += (dim - k - 1) * std::log(1 - pow(partial_correlations[k][i],2) );
        }
    }
    ln_jacobian *= -1/2;
    
    // choose an index
    size_t indexa = size_t( rng->uniform01() * dim );
    size_t indexb = size_t( rng->uniform01() * dim );
    
    // make sure we don't get a diagonal
    while (indexb == indexa)
    {
        indexb = size_t( rng->uniform01() * dim );
    }

    // get the chosen partial correlation
    double current_value = partial_correlations[indexa][indexb];

    // transform the current value from [-1, 1] to [0, 1]
    current_value = (current_value + 1.0) / 2.0;

    // draw new rates and compute the hastings ratio at the same time
    double a = alpha + 1.0;
    double b = (a - 1.0) / current_value - a + 2.0;
    double new_value = RbStatistics::Beta::rv(a, b, *rng);
    
    // set the value (for both sides of the matrix!)
    double new_value_transformed = new_value * 2.0 - 1.0;
    partial_correlations[indexa][indexb] = new_value_transformed;
    partial_correlations[indexb][indexa] = new_value_transformed;
    
    // now transform the partial correlations back into a correlation matrix
    MatrixReal P(dim, dim, 1.0);
    for (int k = 0; k < dim - 1; ++k)
    {
        for (int i = k + 1; i < dim; ++i)
        {
            double p = partial_correlations[k][i]; // initial value for the product-moment correlation
            for (int l = k - 1; l >= 0; --l)
            {
                p = p * pow( ( 1.0 - pow(partial_correlations[l][i], 2) ) * ( 1.0 - pow(partial_correlations[l][k], 2) ), 0.5) + partial_correlations[l][i] * partial_correlations[l][k];
            }
            P[i][k] = p;
            P[k][i] = p;
        }
    }

    // now set the value of the variable
    MatrixReal& current_matrix = variable->getValue();
    for (size_t i = 0; i < (dim - 1); ++i)
    {
        for (size_t j = i+1; j < dim; ++j)
        {
            current_matrix[i][j] = P[i][j];
            current_matrix[j][i] = P[i][j];
            variable->addTouchedElementIndex(i);
            variable->addTouchedElementIndex(j);
        }
    }
    variable->touch(true);

    double ln_hastings_ratio = 0.0;
    try
    {
        // compute the Hastings ratio
        double forward = RbStatistics::Beta::lnPdf(a, b, new_value);
        double new_a = alpha + 1.0;
        double new_b = (a - 1.0) / new_value - a + 2.0;
        double backward = RbStatistics::Beta::lnPdf(new_a, new_b, current_value);
        ln_hastings_ratio = backward - forward;
    }
	catch (RbException &e)
    {
        ln_hastings_ratio = RbConstants::Double::neginf;
    }
    
    return ln_hastings_ratio - ln_jacobian;
    
}


/**
 *
 */
void CorrelationMatrixPartialElementBetaProposal::prepareProposal( void )
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
void CorrelationMatrixPartialElementBetaProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void CorrelationMatrixPartialElementBetaProposal::undoProposal( void )
{
    
    // swap current value and stored value
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
void CorrelationMatrixPartialElementBetaProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    
}


void CorrelationMatrixPartialElementBetaProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void CorrelationMatrixPartialElementBetaProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }
    
}

