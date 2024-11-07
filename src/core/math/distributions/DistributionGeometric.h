/**
 * @file DistributionGeometric.h
 * This file contains the functions of the geometric distribution.
 */


#ifndef DistributionGeometric_H
#define DistributionGeometric_H

#include <cstdint>

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Geometric {
        
            double                      pdf(std::int64_t n, double p);                          /*!< Geometric(p) probability density */
            double                      pdf(std::int64_t n, double p, bool asLog);              /*!< Geometric(p) probability density */
            double                      lnPdf(std::int64_t n, double p);                        /*!< Geometric(p) log_e probability density */
            double                      cdf(std::int64_t n, double p);                          /*!< Geometric(p) cumulative probability */
            std::int64_t                quantile(double q, double p);                           /*!< Geometric(p) quantile */
            std::int64_t                rv(double p, RandomNumberGenerator& rng);               /*!< Geometric(p) random variable */
        }
    }
}

#endif
