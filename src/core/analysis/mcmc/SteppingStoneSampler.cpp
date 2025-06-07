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


SteppingStoneSampler::SteppingStoneSampler(const std::string &fn, const std::string &pn, const std::string &ln,  const std::string &del) : MarginalLikelihoodEstimator(fn, pn, ln, del)
{
    
}



SteppingStoneSampler::~SteppingStoneSampler()
{
    
}



SteppingStoneSampler* SteppingStoneSampler::clone( void ) const
{
    return new SteppingStoneSampler(*this);
}


double SteppingStoneSampler::marginalLikelihood( void ) const
{
    
    double marginal = 0.0;
    
    if ( process_active == true )
    {
        for (size_t i = 1; i < powers.size(); ++i)
        {
            
            size_t samplesPerPath = likelihoodSamples[i].size();
            double max = likelihoodSamples[i][0];
            for (size_t j = 1; j < samplesPerPath; ++j)
            {
                if (max < likelihoodSamples[i][j])
                {
                    max = likelihoodSamples[i][j];
                }
            }
        
            // mean( exp(samples-max)^(beta[k-1]-beta[k]) )     or
            // mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) )
            double mean = 0.0;
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += exp( (likelihoodSamples[i][j] - max)*(powers[i - 1] - powers[i]) ) / samplesPerPath;
            }
        
            marginal += log(mean) + (powers[i - 1] - powers[i])*max;
        
        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return marginal;
}


// The following is adapted from TraceNumeric::update()
double SteppingStoneSampler::getESS(const std::vector<double> values) const
{
    size_t MAX_LAG = 1000;
    size_t samples = values.size();
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
    
    // get the mean
    double mean = 0.0;
    for (size_t i = 0; i < samples; ++i)
    {
        mean += values[i];
    }
    mean /= samples;
    
    double* gammaStat = new double[maxLag];
    for (size_t j = 0; j < maxLag; j++)
    {
        gammaStat[j] = 0;
    }
    double varStat = 0.0;
    
    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values[j] - mean;
            double del2 = values[j + lag] - mean;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            } else {
                maxLag = lag;
            }
        }
    }
    
    // autocorrelation time
    double act = varStat / gammaStat[0];

    // effective sample size
    double ess = samples / act;
    
    return ess;
}


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
