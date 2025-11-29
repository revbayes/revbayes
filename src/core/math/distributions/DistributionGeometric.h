/**
 * @file DistributionGeometric.h
 * This file contains the functions of the geometric distribution.
 */


#ifndef DistributionGeometric_H
#define DistributionGeometric_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Geometric {
        
            double                      pdf(long n, double p);                                   /*!< Geometric(p) probability density */
            double                      pdf(long n, double p, bool asLog);                       /*!< Geometric(p) probability density */
            double                      lnPdf(long n, double p);                                 /*!< Geometric(p) log_e probability density */
            double                      cdf(long n, double p);                                   /*!< Geometric(p) cumulative probability */
            long                         quantile(double q, double p);                           /*!< Geometric(p) quantile */
            long                         rv(double p, RandomNumberGenerator& rng);               /*!< Geometric(p) random variable */
        }
    }
}

#endif
