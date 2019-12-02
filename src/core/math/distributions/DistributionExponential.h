/**
 * @file DistributionExponential
 * This file contains the functions of the exponential distribution.
 *
 * @brief Implementation of the exponential distribution.
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
