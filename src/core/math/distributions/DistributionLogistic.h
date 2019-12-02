/**
 * @file DistributionLogistic
 * This file contains the functions of the logistic distribution.
 *
 * @brief Implementation of the logistic distribution.
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


#ifndef DistributionLogistic_H
#define DistributionLogistic_H

#include <vector>

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace Logistic {
            double      pdf(double x);                                          /*!< Logistic(x) probability density */
            double      pdf(double location,  double scale,  double x);         /*!< Logistic(x) probability density */
            double      lnPdf(double location,  double scale,  double x);       /*!< Logistic(x) log_e probability density */
            double      cdf(double location,  double scale,  double x);         /*!< Logistic(location, scale) cumulative probability */
            double      quantile(double p);                                                     /*!< Logistic quantile */
            double      quantile(double location, double scale, double p);                            /*!< Logistic(location, scale) quantile */
            double      rv(double location, double scale, RandomNumberGenerator& rng);         /*!< Logistic(x) random variable */
        }
    }
}

#endif
