/**
 * @file DistributionChisq
 * This file contains the functions of the chi-square distribution.
 *
 * @brief Implementation of the chi-square distribution.
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


#ifndef DistributionChisq_H
#define DistributionChisq_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace ChiSquare {
        
            double                      pdf(double df, double x);                                              /*!< Chi-Square(v) probability density */
            double                      lnPdf(std::int64_t df, double x);                                              /*!< Chi-Square(v) log_e probability density */
            double                      lnPdf(double df, double x);                                            /*!< Chi-Square(v) log_e probability density */
            double                      cdf(double df, double x);                                              /*!< Chi-Square(v) cumulative probability */
            double                      quantile(double prob, double df);                                      /*!< Chi-Square(v) quantile */
            double                      rv(std::int64_t df, RandomNumberGenerator& rng);                               /*!< Chi-Square(v) random variable */
            double                      rv(double df, RandomNumberGenerator& rng);                             /*!< Chi-Square(v) random variable */
        }
    }
}

#endif
