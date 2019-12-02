/**
 * @file DistributionLaplace
 * This file contains the functions of the normal distribution.
 *
 * @brief Implementation of the normal distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @since 2017-10-19
 *
 * $Id$
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
