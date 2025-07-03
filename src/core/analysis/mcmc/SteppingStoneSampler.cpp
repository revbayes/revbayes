#include "SteppingStoneSampler.h"

#include <cstddef>
#include <cmath>
#include <vector>
#include <iostream>

#include "Cloneable.h"
#include "RlUserInterface.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;


SteppingStoneSampler::SteppingStoneSampler(const std::string &fn, const std::string &pn, const std::string &ln, const std::string &del) : MarginalLikelihoodEstimator(fn, pn, ln, del)
{
    
}



SteppingStoneSampler::~SteppingStoneSampler()
{
    
}



SteppingStoneSampler* SteppingStoneSampler::clone( void ) const
{
    return new SteppingStoneSampler(*this);
}


/**
 * @param power_vec Vector of k powers used in the power posterior analysis that generated the samples
 * @param lnl_vec Nested vector of n log likelihoods sampled by the power posterior analysis for each of the k powers
 * @return Stepping-stone marginal likelihood
 */
double SteppingStoneSampler::marginalLikelihoodGeneral(std::vector<double> power_vec, std::vector< std::vector<double> > lnl_vec) const
{
    
    double marginal = 0.0;
    
    if ( process_active == true )
    {
        for (size_t i = 1; i < power_vec.size(); ++i)
        {
            
            size_t samplesPerPath = lnl_vec[i].size();
            double max = lnl_vec[i][0];
            for (size_t j = 1; j < samplesPerPath; ++j)
            {
                if (max < lnl_vec[i][j])
                {
                    max = lnl_vec[i][j];
                }
            }
        
            // mean( exp(samples-max)^(beta[k-1]-beta[k]) )     or
            // mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) )
            double mean = 0.0;
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += exp( (lnl_vec[i][j] - max)*(power_vec[i - 1] - power_vec[i]) ) / samplesPerPath;
            }
        
            marginal += log(mean) + (power_vec[i - 1] - power_vec[i])*max;
        
        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return marginal;
}


/**
 * Calculate the standard error of the marginal likelihood estimate using delta approximation, as suggested by Xie et al.
 * (2011; Syst. Biol. 60(2): 150--160). The mathematics follows the equations given by Xie et al. (2011: 153), and the specific
 * implementation (including some of the variable names) is inspired by -- and tested against -- mcmc3r::stepping.stones()
 * (see https://github.com/dosreislab/mcmc3r/blob/main/R/marginal-lhd.R), a function written by Mario dos Reis
 * for his companion R package for MCMCTree. Like the latter, this function displays a warning when the delta approximation
 * does not work well. We collect these warnings and present them to the user after the calculation has been completed.
 *
 * @return The approximate standard error
 */
double SteppingStoneSampler::standardError( void ) const
{
    double vmlnl = 0.0;
    
    if ( process_active == true )
    {
        std::stringstream warnings;
        size_t warning_counter = 0;
        
        for (size_t i = 1; i < powers.size(); ++i)
        {
            
            size_t samplesPerPath = likelihoodSamples[i].size();
            
            // compute the maximum of the log likelihood values
            double max = likelihoodSamples[i][0];
            for (size_t j = 1; j < samplesPerPath; ++j)
            {
                if (max < likelihoodSamples[i][j])
                {
                    max = likelihoodSamples[i][j];
                }
            }
            
            /* For numerical stability, we first subtract the maximum sampled log likelihood from the rest.
             * We then multiply the difference by the difference of successive powers, and exponentiate the result.
             * This is equivalent to dividing by the largest sampled likelihood term, and raising the result to
             * the difference of successive powers (Xie et al. 2011: 153).
             */
            std::vector<double> Ls(samplesPerPath);
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                Ls[j] = exp( (likelihoodSamples[i][j] - max)*(powers[i - 1] - powers[i]) );
            }
            
            // get the variance and the mean of the resulting values
            double mean_Ls = 0.0;
            double var_Ls = 0.0;
            
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean_Ls += Ls[j];
            }
            mean_Ls /= samplesPerPath;
            
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                var_Ls += (Ls[j] - mean_Ls)*(Ls[j] - mean_Ls);
            }
            var_Ls /= (int(samplesPerPath) - 1);
            
            // get the effective sample size and divide the variance by it
            double ess = getESS(Ls);
            double vzr = var_Ls / ess;
            
            if (vzr / (mean_Ls*mean_Ls) > 0.1)
            {
                warning_counter++;
                if ( warning_counter == 1 )
                {
                    warnings << "Warning: Standard error approximation unreliable as var(r_k)/r_k^2 > 0.1 for the following powers:\n" << powers[i];
                }
                else
                {
                    warnings << ", " << powers[i];
                }
            }
            
            vmlnl += vzr / (mean_Ls*mean_Ls);
            
        }
        
        if ( not warnings.str().empty() )
        {
            warnings << ".\n";
            RBOUT( warnings.str() );
        }
    }
    
#ifdef RB_MPI
    MPI_Bcast(&vmlnl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return sqrt(vmlnl);
}
