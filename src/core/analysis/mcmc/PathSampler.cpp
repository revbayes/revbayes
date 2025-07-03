#include "PathSampler.h"

#include <cstddef>
#include <cmath>
#include <vector>

#include "Cloneable.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;



PathSampler::PathSampler(const std::string &fn, const std::string &pn, const std::string &ln, const std::string &del) : MarginalLikelihoodEstimator(fn, pn, ln, del)
{
    
}



PathSampler::~PathSampler()
{
    
}



PathSampler* PathSampler::clone( void ) const
{
    return new PathSampler(*this);
}


/**
 * @param power_vec Vector of k powers used in the power posterior analysis that generated the samples
 * @param lnl_vec Nested vector of n log likelihoods sampled by the power posterior analysis for each of the k powers
 * @return Path-sampling marginal likelihood
 */
double PathSampler::marginalLikelihoodGeneral(std::vector<double> power_vec, std::vector< std::vector<double> > lnl_vec) const
{
    
    double marginal = 0.0;
    
    if ( process_active == true )
    {
        // create a vector for the mean log-likelihood values per power posterior
        std::vector<double> pathValues;
    
        // iterate over all powers
        for (size_t i = 0; i < power_vec.size(); ++i)
        {
            // compute the mean for this power
            double mean = 0.0;
        
            // get the number of samples for this power posterior analysis
            size_t samplesPerPath = lnl_vec[i].size();
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += lnl_vec[i][j];
            }
        
            // store the mean
            pathValues.push_back( mean / samplesPerPath );
        
        }
    
        // now we can compute the marginal likelihood
        // the method uses the trapezoidal rule for numerical integration
        for (size_t i = 0; i+1 < pathValues.size(); ++i)
        {
            marginal += (pathValues[i] + pathValues[i + 1])*(power_vec[i] - power_vec[i + 1])/2.0;
        }
        
    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return marginal;
}


double PathSampler::standardError( void ) const
{
    double vmlnl = 0.0;
    
    if ( process_active == true )
    {
        // create vectors for the variance and effective sample size of log-likelihood values per power posterior
        std::vector<double> var_vect( powers.size() );
        std::vector<double> ess_vect( powers.size() );
        
        // iterate over all powers
        for (size_t i = 0; i < powers.size(); ++i)
        {
            size_t samplesPerPath = likelihoodSamples[i].size();
            
            // get the mean and variance for this power
            double mean = 0.0;
            double var = 0.0;
            
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += likelihoodSamples[i][j];
            }
            mean /= samplesPerPath;
            
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                var += (likelihoodSamples[i][j] - mean)*(likelihoodSamples[i][j] - mean);
            }
            var /= (int(samplesPerPath) - 1);
            
            var_vect[i] = var;
            
            // get the effective sample size
            ess_vect[i] = getESS( likelihoodSamples[i] );
        }
        
        // for the first sampling variance term (indexed 0), the weight is computed as (beta_1 - beta_0)/2
        vmlnl += (var_vect[0] / ess_vect[0])*( (powers[1] - powers[0])*(powers[1] - powers[0]) / 4 );
        
        // for all sampling variance terms except the first and the last one, the weight is computed as
        // (beta_(k + 1) - beta_(k - 1))/2
        size_t K = powers.size() - 1;
        for (size_t i = 1; i < K; ++i)
        {
            double vv = (var_vect[i] / ess_vect[i]);
            double squared_weight = (powers[i + 1] - powers[i - 1])*(powers[i + 1] - powers[i - 1]) / 4;
            vmlnl += vv*squared_weight;
        }
        
        // for the last sampling variance term (indexed K), the weight is computed as (beta_K - beta_(K - 1))/2
        vmlnl += (var_vect[K] / ess_vect[K])*( (powers[K] - powers[K - 1])*(powers[K] - powers[K - 1]) / 4 );
    }
    
#ifdef RB_MPI
    MPI_Bcast(&vmlnl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return sqrt(vmlnl);
}
