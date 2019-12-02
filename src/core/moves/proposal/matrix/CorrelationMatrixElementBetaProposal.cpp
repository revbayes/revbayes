#include "CorrelationMatrixElementBetaProposal.h"
#include "DistributionBeta.h"
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
CorrelationMatrixElementBetaProposal::CorrelationMatrixElementBetaProposal( StochasticNode<MatrixReal> *n, double a, double p) : Proposal(p),
    variable( n ),
    alpha( a ),
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
void CorrelationMatrixElementBetaProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixElementBetaProposal* CorrelationMatrixElementBetaProposal::clone( void ) const
{
    
    return new CorrelationMatrixElementBetaProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixElementBetaProposal::getProposalName( void ) const
{
    static std::string name = "CorrelationMatrixSingleElementBeta";
    
    return name;
}


double CorrelationMatrixElementBetaProposal::getProposalTuningParameter( void ) const
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
double CorrelationMatrixElementBetaProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MatrixReal& v = variable->getValue();
    // choose an index
    indexa = size_t( rng->uniform01() * v.getNumberOfRows() );
    indexb = size_t( rng->uniform01() * v.getNumberOfColumns() );
    
    // make sure we don't get a diagonal
    while (indexb == indexa)
    {
        indexb = size_t( rng->uniform01() * v.getNumberOfColumns() );
    }
    
    // copy the current value
    storedValue = v[indexa][indexb];
    double current_value = storedValue;
    
    // transform the current value from [-1, 1] to [0, 1]
    current_value = (current_value + 1.0) / 2.0;
    
    // draw new rates and compute the hastings ratio at the same time
    double a = alpha + 1.0;
    double b = (a - 1.0) / current_value - a + 2.0;
//    double a = alpha * current_value;
//    double b = alpha * (1 - current_value);
    double new_value = RbStatistics::Beta::rv(a, b, *rng);

    // set the value (for both sides of the matrix!)
    double new_value_transformed = new_value * 2.0 - 1.0;
    v[indexa][indexb] = new_value_transformed;
    v[indexb][indexa] = new_value_transformed;
    
    variable->addTouchedElementIndex(indexa);
    variable->addTouchedElementIndex(indexb);
    
    double ln_Hastings_ratio = 0.0;
    try
    {
        // compute the Hastings ratio
        double forward = RbStatistics::Beta::lnPdf(a, b, new_value);
        double new_a = alpha + 1.0;
        double new_b = (a - 1.0) / new_value - a + 2.0;
//        double new_a = alpha * new_value;
//        double new_b = alpha * (1 - new_value);
        double backward = RbStatistics::Beta::lnPdf(new_a, new_b, current_value);
        ln_Hastings_ratio = backward - forward;
    }
	catch (RbException &e)
    {
        ln_Hastings_ratio = RbConstants::Double::neginf;
    }
    
    return ln_Hastings_ratio;
    
}


/**
 *
 */
void CorrelationMatrixElementBetaProposal::prepareProposal( void )
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
void CorrelationMatrixElementBetaProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void CorrelationMatrixElementBetaProposal::undoProposal( void )
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
void CorrelationMatrixElementBetaProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    
}


void CorrelationMatrixElementBetaProposal::setProposalTuningParameter(double tp)
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
void CorrelationMatrixElementBetaProposal::tune( double rate )
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

