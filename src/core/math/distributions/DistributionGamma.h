/**
 * @file DistributionGamma.h
 * This file contains the functions of the gamma distribution.
 *
 */


#ifndef DistributionGamma_H
#define DistributionGamma_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Gamma {
        
            double                      pdf(double a, double b, double x);                                    /*!< Gamma(a,b) probability density */
            double                      pdf(double a, double b, double x, bool isLog);                        /*!< Gamma(a,b) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< Log of the Gamma(a,b) probability density */  
            double                      cdf(double a, double b, double x);                                    /*!< Gamma(a,b) cumulative probability */
            double                      quantile(double a, double b, double p);                               /*!< Gamma(a,b) quantile */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< Gamma(a,b) random variable */
        
        }
    }
}

#endif
