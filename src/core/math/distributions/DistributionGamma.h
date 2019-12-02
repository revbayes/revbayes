/**
 * @file DistributionGamma
 * This file contains the functions of the gamma distribution.
 *
 * @brief Implementation of the gamma distribution.
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
