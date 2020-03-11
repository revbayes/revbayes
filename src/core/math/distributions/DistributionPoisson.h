/**
 * @file DistributionPoisson.h
 * This file contains the functions of the poisson distribution.
 *
 */


#ifndef DistributionPoisson_H
#define DistributionPoisson_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace Poisson {
        
            double                      pdf(double lambda, int x);                                            /*!< Poisson(lambda) probability */
            double                      lnPdf(double lambda, int x);                                          /*!< Log of the Poisson(lambda) probability */
            double                      cdf(double lambda, int x);                                            /*!< Poisson(lambda) cumulative probability */
            double                      quantile(double lambda, double p);                                    /*!< Poisson(lambda) quantile */
            int                         rv(double lambda, RandomNumberGenerator& rng);                        /*!< Poisson(lambda) random variable */
        }
    }
}

#endif
