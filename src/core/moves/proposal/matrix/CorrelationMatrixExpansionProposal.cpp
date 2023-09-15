#include "CorrelationMatrixExpansionProposal.h"

#include <cstddef>
#include <cmath>
#include <vector>

#include "DistributionInverseWishart.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "Cloneable.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class ContinuousCharacterData; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
CorrelationMatrixExpansionProposal::CorrelationMatrixExpansionProposal( StochasticNode<MatrixReal> *cm, TypedDagNode<RbVector<double> > *v, StochasticNode<ContinuousCharacterData>* mv) : Proposal(),
    correlation_matrix( cm ),
    variance( v ),
    mvbm( mv ),
    stored_matrix( cm->getValue().getNumberOfRows() )
{
    
    PhyloMultivariateBrownianProcessREML* dist = dynamic_cast<PhyloMultivariateBrownianProcessREML *>( &mvbm->getDistribution() );
    if (dist == NULL)
    {
        throw(RbException("CorrelationMatrixExpansion move only works for multivariate phylogenetic Brownian motion models."));
    }
    
    addNode( correlation_matrix );
    addNode( variance );
    addNode( mvbm );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CorrelationMatrixExpansionProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixExpansionProposal* CorrelationMatrixExpansionProposal::clone( void ) const
{
    
    return new CorrelationMatrixExpansionProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixExpansionProposal::getProposalName( void ) const
{
    static std::string name = "MatrixRealSingleElementBetaMove";
    
    return name;
}


double CorrelationMatrixExpansionProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
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
double CorrelationMatrixExpansionProposal::doProposal( void )
{

    stored_matrix = correlation_matrix->getValue();
    stored_matrix.setCholesky(true);
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the standardized contrasts
    PhyloMultivariateBrownianProcessREML* dist         = dynamic_cast<PhyloMultivariateBrownianProcessREML *>( &mvbm->getDistribution() );
    std::vector<std::vector<double> >     contrasts    = dist->getContrasts();
    std::vector<double>                   contrast_sds = dist->getContrastStDevs();
    
    // compute some helpers
    size_t num_contrasts = contrasts.size();
    size_t num_traits    = correlation_matrix->getValue().getNumberOfRows();
    
    // compute the sum of squared values for each trait across contrasts
    std::vector<double> squared_trait_sums(num_traits, 0.0);
    for (size_t i = 0; i < num_contrasts; ++i)
    {
        for (size_t j = 0; j < num_traits; ++j)
        {
            squared_trait_sums[j] += pow(contrasts[i][j], 2.0);
        }
    }
    
    // rescale the contrasts so that the sum of squared contrasts for each trait is 1
    std::vector<double> scale_factor(num_traits, 0.0);
    for (size_t i = 0; i < num_traits; ++i)
    {
        scale_factor[i] = pow(1.0 / squared_trait_sums[i], 0.5);
    }
    
    std::vector<std::vector<double> > transformed_contrasts = contrasts;
    for (size_t i = 0; i < num_contrasts; ++i)
    {
        for (size_t j = 0; j < num_traits; ++j)
        {
            transformed_contrasts[i][j] *= scale_factor[j];
        }
    }
    
    // compute the parameters of the inverse Wishart distribution
    int nu = int(num_contrasts);
    MatrixReal E(num_traits, num_traits);
    MatrixReal V(num_traits, num_traits);
    MatrixReal S(num_traits, num_traits);
    std::vector<double> var = variance->getValue();
    
    for (size_t i = 0; i < num_contrasts; ++i)
    {
        std::vector<double> &these_contrasts  = transformed_contrasts[i];
        double              &this_contrast_sd = contrast_sds[i];
        for (size_t j = 0; j < num_traits; j++)
        {
            for (size_t k = j; k < num_traits; k++)
            {
                E[j][k] = these_contrasts[j] * these_contrasts[k];
                E[k][j] = these_contrasts[j] * these_contrasts[k];
            }
            V[j][j] = 1.0 / (pow(var[j], 0.5) * this_contrast_sd);
        }
        S += V * E * V;
    }
    
    // sample a new variance-covariance matrix from the inverse Wishart distribution
//    MatrixReal new_vcv = RbStatistics::InverseWishart::rvCovariance(S, nu, *rng);
    MatrixReal new_vcv = RbStatistics::InverseWishart::rv(S, nu, *rng);
    
    // compute the scaling matrix, D, necessary to get back to the correlation matrix
    MatrixReal D_inv(num_traits, num_traits);
    for (size_t i = 0; i < num_traits; ++i)
    {
        D_inv[i][i] = 1.0 / pow(new_vcv[i][i], 0.5);
    }
    MatrixReal new_correlation_matrix = D_inv * new_vcv * D_inv;
    new_correlation_matrix.setCholesky(true);
    
    // compute the hastings ratio
    double ln_hastings_ratio = ((num_traits + 1.0) / 2.0) * (new_correlation_matrix.getLogDet() - stored_matrix.getLogDet());
    
    // update the parameters and compute the prior ratio
//    double ln_likelihood_before = mvbm->getLnProbability();
    MatrixReal& current_matrix = correlation_matrix->getValue();
    for (size_t i = 0; i < num_traits; ++i)
    {
        for (size_t j = i; j < num_traits; ++j){
            current_matrix[i][j] = new_correlation_matrix[i][j];
            current_matrix[j][i] = new_correlation_matrix[i][j];
        }
    }
    correlation_matrix->touch(true);
    
    // we need to subtract the likelihood ratio
    double ln_likelihood_ratio = mvbm->getLnProbabilityRatio();
//    double ln_likelihood_after = mvbm->getLnProbability();
    
    // subtract the likelihood ratio because this is supposed to be a Gibbs move
    ln_hastings_ratio -= ln_likelihood_ratio;

//    std::cout << ln_likelihood_before << " -- " << ln_likelihood_after << std::endl;
    
    return ln_hastings_ratio;
    
}


/**
 *
 */
void CorrelationMatrixExpansionProposal::prepareProposal( void )
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
void CorrelationMatrixExpansionProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    ;
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void CorrelationMatrixExpansionProposal::undoProposal( void )
{
    
    MatrixReal& current_matrix = correlation_matrix->getValue();
    size_t size = current_matrix.getNumberOfRows();
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = i; j < size; ++j){
            current_matrix[i][j] = stored_matrix[i][j];
            current_matrix[j][i] = stored_matrix[i][j];
        }
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void CorrelationMatrixExpansionProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == correlation_matrix)
    {
        correlation_matrix = static_cast<StochasticNode<MatrixReal>* >(newN) ;
    }
    else if (oldN == variance)
    {
        variance = static_cast<TypedDagNode<RbVector<double> >* >(newN) ;
    }
    else if (oldN == mvbm)
    {
        mvbm = static_cast<StochasticNode<ContinuousCharacterData>* >(newN) ;
    }
    
}


void CorrelationMatrixExpansionProposal::setProposalTuningParameter(double tp)
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
void CorrelationMatrixExpansionProposal::tune( double rate )
{
    
//    double p = this->targetAcceptanceRate;
//    
//    if ( rate > p )
//    {
//        alpha /= (1.0 + ((rate-p)/(1.0 - p)) );
//    }
//    else
//    {
//        alpha *= (2.0 - rate/p);
//    }
    
}

