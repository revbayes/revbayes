/**
 * @file DistributionKumaraswamy
 * This file contains the functions of the Kumaraswamy distribution.
 *
 * @brief Implementation of the Kumaraswamy distribution.
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


#ifndef DistributionKumaraswamy_H
#define DistributionKumaraswamy_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace Kumaraswamy {

            double                      pdf(double a, double b, double x);                                    /*!< Kumaraswamy(a,b) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< Kumaraswamy(a,b) log_e probability density */
            double                      cdf(double a, double b, double x);                                    /*!< Kumaraswamy(a,b) cumulative probability */
            double                      quantile(double a, double b, double p);                               /*!< Kumaraswamy(a,b) quantile */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< Kumaraswamy(a,b) random variable */
        }
    }
}

#endif
