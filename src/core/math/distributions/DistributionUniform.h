/**
 * @file DistributionUniform
 * This file contains the functions of the uniform distribution.
 *
 * @brief Implementation of the uniform distribution.
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


#ifndef DistributionUniform_H
#define DistributionUniform_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Uniform {
        
            double                      pdf(double x);                                                        /*!< Uniform(0,1) probability density */
            double                      pdf(double a, double b, double x);                                    /*!< Uniform(a,b) probability density */
            double                      lnPdf(double x);                                                      /*!< Log of the Uniform(0,1) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< Log of the Uniform(a,b) probability density */
            double                      cdf(double x);                                                        /*!< Uniform(0,1) cumulative probability */
            double                      cdf(double a, double b, double x);                                    /*!< Uniform(a,b) cumulative probability */
            inline                      double quantile(double p);                                            /*!< Uniform(0,1) quantile */
            double                      quantile(double a, double b, double p);                               /*!< Uniform(a,b) quantile */
            double                      rv(RandomNumberGenerator& rng);                                       /*!< Uniform(0,1) random variable */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< Uniform(a,b) random variable */
        }
    }
}

#endif
