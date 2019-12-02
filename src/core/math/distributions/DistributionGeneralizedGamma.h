/**
 * @file DistributionGamma
 * This file contains the functions of the generalized gamma distribution.
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


#ifndef DistributionGeneralizedGamma_H
#define DistributionGeneralizedGamma_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace GeneralizedGamma {

            double                      pdf(double a, double b, double l, double x);                                    /*!< Gamma(a,b) probability density */
            double                      lnPdf(double a, double b, double l, double x);                                  /*!< Log of the Gamma(a,b) probability density */
            double                      cdf(double a, double b, double l, double x);                                    /*!< Gamma(a,b) cumulative probability */
            double                      quantile(double a, double b, double l, double p);                               /*!< Gamma(a,b) quantile */
            double                      rv(double a, double b, double l, RandomNumberGenerator& rng);                   /*!< Gamma(a,b) random variable */

        }
    }
}

#endif
