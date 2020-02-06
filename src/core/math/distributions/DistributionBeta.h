/**
 * @file DistributionBeta.h
 *
 * @brief Implementation of the beta distribution.
 *
 */


#ifndef DistributionBeta_H
#define DistributionBeta_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    namespace RbStatistics {
    
        namespace Beta {
        
            double                      pdf(double a, double b, double x);                                    /*!< Beta(a,b) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< Beta(a,b) log_e probability density */
            double                      cdf(double a, double b, double x);                                    /*!< Beta(a,b) cumulative probability */
            double                      quantile(double a, double b, double p);                               /*!< Beta(a,b) quantile */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< Beta(a,b) random variable */
        }
    }
}

#endif
