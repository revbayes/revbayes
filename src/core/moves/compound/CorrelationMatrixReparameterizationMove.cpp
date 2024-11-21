#include "CorrelationMatrixReparameterizationMove.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "DistributionInverseWishart.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class ContinuousCharacterData; }

using namespace RevBayesCore;

CorrelationMatrixReparameterizationMove::CorrelationMatrixReparameterizationMove(StochasticNode<MatrixReal> *cm, TypedDagNode<RbVector<double> > *v, StochasticNode<ContinuousCharacterData>* mv, bool t, double w) : AbstractMove( w, t),
    correlation_matrix( cm ),
    variance( v ),
    mvbm( mv ),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 ),
    stored_matrix( cm->getValue().getNumberOfRows() )
{
    
    PhyloMultivariateBrownianProcessREML* dist = dynamic_cast<PhyloMultivariateBrownianProcessREML *>( &mvbm->getDistribution() );
    if (dist == NULL)
    {
        throw(RbException("CorrelationMatrixReparameterization move only works for multivariate phylogenetic Brownian motion models."));
    }
    
    addNode( correlation_matrix );
    addNode( variance );
    addNode( mvbm );
    
}


/**
 * Basic destructor doing nothing.
 */
CorrelationMatrixReparameterizationMove::~CorrelationMatrixReparameterizationMove( void )
{
    // nothing special to do
    // everything should be taken care of in the base class
    
}




/* Clone object */
CorrelationMatrixReparameterizationMove* CorrelationMatrixReparameterizationMove::clone( void ) const
{
    
    return new CorrelationMatrixReparameterizationMove( *this );
}



const std::string& CorrelationMatrixReparameterizationMove::getMoveName( void ) const
{
    
    static std::string name = "CorrelationMatrixReparameterizationMove";
    
    return name;
}


double CorrelationMatrixReparameterizationMove::getMoveTuningParameter( void ) const
{
    // this move has no tuning parameter
    return RbConstants::Double::nan;
}


size_t CorrelationMatrixReparameterizationMove::getNumberAcceptedCurrentPeriod( void ) const
{
    return num_accepted_current_period;
}


size_t CorrelationMatrixReparameterizationMove::getNumberAcceptedTotal( void ) const
{
    return num_accepted_total;
}


/** Perform the move */
void CorrelationMatrixReparameterizationMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
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
        S += (V * E * V);
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
    double old_ln_lik = dist->computeLnProbability();
    
    double ln_hastings_ratio = ((num_traits + 1.0) / 2.0) * (new_correlation_matrix.getLogDet() - stored_matrix.getLogDet());
    
    // update the parameters and compute the prior ratio
    correlation_matrix->setValue( new MatrixReal(new_correlation_matrix) );
    double ln_prior_ratio = correlation_matrix->getLnProbabilityRatio();
    
    double new_ln_lik = dist->computeLnProbability();

    // compute the acceptance ratio
    double ln_acceptance_ratio = pHeat * prHeat * ln_prior_ratio + ln_hastings_ratio;
//    ln_acceptance_ratio = -500.0;
    double acceptance_ratio = exp(ln_acceptance_ratio);
    

    if (ln_acceptance_ratio > 0.0)
    {
        num_accepted_total++;
        num_accepted_current_period++;
        
        correlation_matrix->keep();
        mvbm->keep();
        double lnp = dist->computeLnProbability();
        std::cout << "reject " << acceptance_ratio << " -- " << old_ln_lik << " -- " << new_ln_lik << " -- " << lnp << std::endl;
    }
    else if (ln_acceptance_ratio < -300.0)
    {
        reject();
//        correlation_matrix->restore();
//        mvbm->restore();
        double lnp = dist->computeLnProbability();
        std::cout << "accept " << acceptance_ratio << " -- " << old_ln_lik << " -- " << new_ln_lik << " -- " << lnp << std::endl;
    }
    else
    {
        double r = exp(ln_acceptance_ratio);
        double u = GLOBAL_RNG->uniform01();
        if (u < r)
        {
            num_accepted_total++;
            num_accepted_current_period++;
            
            correlation_matrix->keep();
            mvbm->keep();
            double lnp = dist->computeLnProbability();
            std::cout << "reject " << acceptance_ratio << " -- " << old_ln_lik << " -- " << new_ln_lik << " -- " << lnp << std::endl;
        }
        else
        {
            reject();
            //        correlation_matrix->restore();
            //        mvbm->restore();
            double lnp = dist->computeLnProbability();
            std::cout << "accept " << acceptance_ratio << " -- " << old_ln_lik << " -- " << new_ln_lik << " -- " << lnp << std::endl;
        }
    }

}


void CorrelationMatrixReparameterizationMove::printSummary(std::ostream &o, bool current_period) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the name
    const std::string &n = getMoveName();
    size_t spaces = 40 - (n.length() > 40 ? 40 : n.length());
    o << n;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";
    
    // print the DagNode name
    const std::string &dn_name = (*nodes.begin())->getName();
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";
    
    // print the weight
    int w_length = 4;
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i)
    {
        o << " ";
    }
    o << weight;
    o << " ";
    
    size_t num_tried = num_tried_total;
    size_t num_accepted = num_accepted_total;
    if (current_period == true)
    {
        num_tried = num_tried_current_period;
        num_accepted = num_accepted_current_period;
    }
    
    // print the number of tries
    int t_length = 9;
    if (num_tried > 0) t_length -= (int)log10(num_tried);
    for (int i = 0; i < t_length; ++i)
    {
        o << " ";
    }
    o << num_tried;
    o << " ";
    
    // print the number of accepted
    int a_length = 9;
    if (num_accepted > 0) a_length -= (int)log10(num_accepted);
    
    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_accepted;
    o << " ";
    
    // print the acceptance ratio
    double ratio = num_accepted / (double)num_tried;
    if (num_tried == 0) ratio = 0;
    int r_length = 5;
    
    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";
    
//    proposal->printParameterSummary( o );
    
    o << std::endl;
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
}


void CorrelationMatrixReparameterizationMove::reject( void )
{
    correlation_matrix->setValue( new MatrixReal(stored_matrix) );
//    correlation_matrix->setValue( stored_matrix.clone() );
}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void CorrelationMatrixReparameterizationMove::resetMoveCounters( void )
{
    num_accepted_current_period = 0;
}


void CorrelationMatrixReparameterizationMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


void CorrelationMatrixReparameterizationMove::setMoveTuningParameter(double tp)
{
    // this move has no tuning parameter: nothing to do
}


void CorrelationMatrixReparameterizationMove::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_accepted_current_period = na;
}


void CorrelationMatrixReparameterizationMove::setNumberAcceptedTotal( size_t na )
{
    num_accepted_total = na;
}


void CorrelationMatrixReparameterizationMove::tune( void )
{

//    double rate = numAccepted / double(num_tried);
//    
//    if ( rate > 0.44 )
//    {
//        delta /= (1.0 + ((rate-0.44)/0.56) );
//    }
//    else
//    {
//        delta *= (2.0 - rate/0.44 );
//    }
    
}


