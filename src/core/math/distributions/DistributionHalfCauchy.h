/**
 * @file DistributionChauchy
 * This file contains the functions of the chauchy distribution.
 *
 * @brief Implementation of the half-cauchy distribution (C+).
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Andrew Magee
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#ifndef DistributionHalfCauchy_H
#define DistributionHalfCauchy_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace HalfCauchy {
            double                      pdf(double x);                                                    /*!< Cauchy(a,b) probability density */
            double                      pdf(double location, double scale, double x);                     /*!< Cauchy(a,b) probability density */
            double                      lnPdf(double x);                                                  /*!< Cauchy(a,b) log_e probability density */
            double                      lnPdf(double location, double scale, double x);                   /*!< Cauchy(a,b) log_e probability density */
            double                      cdf(double x);                                                    /*!< Cauchy(a,b) cumulative probability */
            double                      cdf(double location, double scale, double x);                     /*!< Cauchy(a,b) cumulative probability */
            double                      quantile(double p);                                               /*!< Cauchy(a,b) quantile */
            double                      quantile(double location, double scale, double p);                /*!< Cauchy(a,b) quantile */
            double                      rv(RandomNumberGenerator& rng);                                   /*!< Cauchy(a,b) random variable */
            double                      rv(double location, double scale, RandomNumberGenerator& rng);    /*!< Cauchy(a,b) random variable */
	
        }
    }
}

#endif
