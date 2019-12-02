/**
 * @file DistributionNormal
 * This file contains the functions of the normal distribution.
 *
 * @brief Implementation of the normal distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#ifndef DistributionNormal_H
#define DistributionNormal_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Normal {
        
            double                      pdf(double x);                                                          /*!< Normal(0,1) probability density */
            double                      pdf(double mu, double sigma, double x);                                 /*!< Normal(mu,sigma) probability density */
            double                      pdf(double mu, double sigma, double x, double min, double max);         /*!< Truncated Normal probability density */
            double                      lnPdf(double x);                                                        /*!< Log of the Normal(0,1) probability density */
            double                      lnPdf(double mu, double sigma, double x);                               /*!< Log of the Normal(mu,sigma) probability density */
            double                      lnPdf(double mu, double sigma, double x, double min, double max);       /*!< Log of the truncated normal probability density */
            double                      cdf(double x);                                                          /*!< Normal(0,1) cumulative probability */
            double                      cdf(double mu, double sigma, double x);                                 /*!< Normal(mu,sigma) cumulative probability */
            double                      cdf(double mu, double sigma, double x, double min, double max);         /*!< Truncated normal cumulative probability */
            double                      quantile(double p);                                                     /*!< Normal(0,1) quantile */
            double                      quantile(double mu, double sigma, double p);                            /*!< Normal(mu,sigma) quantile */
            double                      quantile(double mu, double sigma, double p, double min, double max);    /*!< Truncated normal quantile */
            double                      rv(RandomNumberGenerator& rng);                                         /*!< Normal(0,1) random variable */
            double                      rv(double mu, double sigma, RandomNumberGenerator& rng);                /*!< Normal(mu,sigma) random variable */
            double                      rv(double mu, double sigma, double min, double max, RandomNumberGenerator& rng); /*!< Truncated normal random variable */
        }
    }
}

#endif
