/**
 * @file DistributionLognormal
 * This file contains the functions of the lognormal distribution.
 *
 * @brief Implementation of the lognormal distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 * Created by sibonli on 8/17/11.
 */


#ifndef DistributionLognormal_H
#define DistributionLognormal_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Lognormal {
        
            double                      pdf(double x);                                                          /*!< Lognormal(0,1) probability density */
            double                      pdf(double m, double s, double x);                                 /*!< Lognormal(m,s) probability density */
            double                      lnPdf(double x);                                                        /*!< Log of the Lognormal(0,1) probability density */
            double                      lnPdf(double m, double s, double x);                               /*!< Log of the Lognormal(m,s) probability density */
            double                      cdf(double x);                                                          /*!< Lognormal(0,1) cumulative probability */
            double                      cdf(double m, double s, double x);                                 /*!< Lognormal(m,s) cumulative probability */
            double                      quantile(double p);                                                     /*!< Lognormal(0,1) quantile */
            double                      quantile(double m, double s, double p);                            /*!< Lognormal(m,s) quantile */
            double                      rv(RandomNumberGenerator& rng);                                         /*!< Lognormal(0,1) random variable */
            double                      rv(double m, double s, RandomNumberGenerator& rng);                /*!< Lognormal(m,s) random variable */
        }
    }
}

#endif
