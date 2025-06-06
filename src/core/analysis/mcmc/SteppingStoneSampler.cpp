#include "SteppingStoneSampler.h"

#include <cstddef>
#include <cmath>
#include <vector>

#include "Cloneable.h"

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
                mean += exp( (likelihoodSamples[i][j]-max)*(powers[i-1]-powers[i]) ) / samplesPerPath;
            }
        
            marginal += log(mean) + (powers[i-1]-powers[i])*max;
        
        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return marginal;
}


double SteppingStoneSampler::standardError( void ) const
{
    size_t MAX_LAG = 1000;
    
    if ( process_active == true )
    {
        double vzr = 0.0;
        double zr = 0.0;
        
        for (size_t i = 1; i < powers.size(); ++i)
        {
            
            size_t samplesPerPath = likelihoodSamples[i].size();
            
            // compute the mean and maximum of the log likelihood values
            double mean = 0.0;
            double max = likelihoodSamples[i][0];
            
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += likelihoodSamples[i][j];
            }
            mean /= samplesPerPath;
            
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
                Ls[i] = exp(  );
            
            // The following is adapted from TraceNumeric::update()
            size_t maxLag = (samplesPerPath - 1 < MAX_LAG ? samplesPerPath - 1 : MAX_LAG);
            
            double* gammaStat = new double[maxLag];
            for (size_t i = 0; i < maxLag; i++)
            {
                gammaStat[i] = 0;
            }
            double varStat = 0.0;
            
            for (size_t lag = 0; lag < maxLag; lag++) {
                for (size_t j = 0; j < samplesPerPath - lag; j++) {
                    double del1 = likelihoodSamples[i][j] - mean;
                    double del2 = likelihoodSamples[i][j + lag] - mean;
                    gammaStat[lag] += (del1 * del2);
                }

                gammaStat[lag] /= ((double) (samplesPerPath - lag));

                if (lag == 0) {
                    varStat = gammaStat[0];
                } else if (lag % 2 == 0) {
                    if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                        varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
                    }
                    else {
                        maxLag = lag;
                    }
                }
            }
            
            // auto correlation time
            double act = varStat / gammaStat[0];

            // effective sample size
            double ess = samplesPerPath / act;
            
            vzr += var / ess;
            
        }
    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return sqrt(vmlnl);
}
