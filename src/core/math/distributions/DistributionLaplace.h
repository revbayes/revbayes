/**
 * @file DistributionLaplace
 * This file contains the functions for the Laplace (or double exponential) distribution.
 *
 * @brief Implementation of the Laplace distribution.
 *
 * The Laplace probability distribution, also referred to as the "Double Exponential"
 * distribution, has probability density
 * @f[ f(x|\mu, b) = {1 \over 2b} \exp \left( - { |x-\mu| \over b} \right) @f]
 * where @f$ \mu @f$ is the location parameter and @f$ b @f$ is a scale parameter.
 * The Laplace distribution is essentially two exponential distributions glued together back-to-back
 * (hence the name, double exponential). The standard Laplace distribution has
 * @f$ \mu=0 @f$ and @f$ b=1 @f$ (i.e., two exponentials with parameter one glued
 * back-to-back and starting at zero).
 *
 */


#ifndef DistributionLaplace_H
#define DistributionLaplace_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Laplace {
        
            double                      pdf(double x);                                                          /*!< Laplace(0,1) probability density */
            double                      pdf(double mu, double scale, double x);                                 /*!< Laplace(mu,scale) probability density */
            double                      lnPdf(double x);                                                        /*!< Log of the Laplace(0,1) probability density */
            double                      lnPdf(double mu, double scale, double x);                               /*!< Log of the Laplace(mu,scale) probability density */
            double                      cdf(double x);                                                          /*!< Laplace(0,1) cumulative probability */
            double                      cdf(double mu, double scale, double x);                                 /*!< Laplace(mu,scale) cumulative probability */
            double                      quantile(double p);                                                     /*!< Laplace(0,1) quantile */
            double                      quantile(double mu, double scale, double p);                            /*!< Laplace(mu,scale) quantile */
            double                      rv(RandomNumberGenerator& rng);                                         /*!< Laplace(0,1) random variable */
            double                      rv(double mu, double scale, RandomNumberGenerator& rng);                /*!< Laplace(mu,scale) random variable */
        }
    }
}

#endif
