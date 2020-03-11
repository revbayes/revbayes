/**
 * @file DistributionExponential.h
 * This file contains the functions of the exponential distribution.
 *
 * $Id$
 */


#ifndef DistributionExponential_H
#define DistributionExponential_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Exponential {
        
            double                      pdf(double lambda, double x);                                         /*!< Exponential(lambda) probability density */
            double                      lnPdf(double lambda, double x);                                       /*!< Exponential(lambda) log_e probability density */
            double                      cdf(double lambda, double x);                                         /*!< Exponential(lambda) cumulative probability */
            double                      quantile(double lambda, double p);                                    /*!< Exponential(lambda) quantile */
            double                      rv(double lambda, RandomNumberGenerator& rng);                        /*!< Exponential(lambda) random variable */
        }
    }
}

#endif
