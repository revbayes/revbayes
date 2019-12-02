#include "PathSampler.h"

#include <stddef.h>
#include <vector>

#include "Cloneable.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;



PathSampler::PathSampler(const std::string &fn, const std::string &pn, const std::string &ln,  const std::string &del) : MarginalLikelihoodEstimator(fn, pn, ln, del)
{
    
}



PathSampler::~PathSampler()
{
    
}



PathSampler* PathSampler::clone( void ) const
{
    return new PathSampler(*this);
}


double PathSampler::marginalLikelihood( void ) const
{
    
    double marginal = 0.0;
    if ( process_active == true )
    {
        // create a vector for the mean log-likelihood values per power posterior
        std::vector<double> pathValues;
    
        // iterate over all powers
        for (size_t i = 0; i < powers.size(); ++i)
        {
            // compute the mean for this power
            double mean = 0.0;
        
            // get the number of samples for this power posterior analysis
            size_t samplesPerPath = likelihoodSamples[i].size();
            for (size_t j = 0; j < samplesPerPath; ++j)
            {
                mean += likelihoodSamples[i][j];
            }
        
            // store the mean
            pathValues.push_back( mean  / samplesPerPath );
        
        }
    
//        // now we can compute the marginal likelihood
//        // the method uses the trapezoidal rule for numerical integration
//        for (size_t i = 0; i < pathValues.size()-1; ++i)
//        {
//            marginal += (pathValues[i] + pathValues[i+1])*(powers[i]-powers[i+1])/2.0;
//        }
     
        // now we compute the marginal likelihood
        // this method uses the generalized simpson's rule (with unequal intervals),
        // equation 11 from:
        // @article{mcnamee1981program,
        //     title={A program to integrate a function tabulated at unequal intervals},
        //     author={McNAMEE, JOHN MICHAEL},
        //     journal={International Journal for Numerical Methods in Engineering},
        //     volume={17},
        //     number={2},
        //     pages={271--279},
        //     year={1981},
        //     publisher={Wiley Online Library}
        // }
        
        // in the case of an odd number of stones, we use the trapezoid rule for
        // one of the integrals. we move backwards from beta=0 to beta=1, so that
        // the trapezoid rule is applied to the stone at beta=1, which we expect
        // to be the "flattest" interval.
        
        int num_powers = (int)pathValues.size();
        for (int i = num_powers - 2; i >= 0; i -= 2)
        {
            
            if (i == 0) {
                
                // if this is the last stone, the just use
                // the trapezoid rule
                
                marginal += (pathValues[i] + pathValues[i+1]) * (powers[i]-powers[i+1]) / 2.0;
                
            } else {
                
                // otherwise, use the generalized simpson's rule
                double x0 = powers[i+1];
                double x1 = powers[i];
                double x2 = powers[i-1];
                
                double y0 = pathValues[i+1];
                double y1 = pathValues[i];
                double y2 = pathValues[i-1];
                
                marginal += (x2 - x0) * (y0 + ( (x2 - x0) / (x1 - x0) ) * (y1 - y0) / 2 );
                marginal += (2.0 * x2 * x2 - x0 * x2 - x0 * x0 + 3 * x0 * x1 - 3 * x1 * x2 ) * ( ((y2 - y1) / (x2 - x1)) - ((y1 - y0) / (x1 - x0)) ) / 6.0 ;
                
            }
            
        }
        
        
    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return marginal;
}

