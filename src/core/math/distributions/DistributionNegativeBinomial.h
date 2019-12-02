/**
 * @file DistributionNegativeBinomial
 * This file contains the functions of the Negative Binomial distribution.
 *
 * @brief Implementation of the Negative Binomial distribution.
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


#ifndef DistributionNegativeBinomial_H
#define DistributionNegativeBinomial_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace NegativeBinomial {
        
            double                      pdf(double r, double p, double x);                                    /*!< NegativeBinomial(r,p) probability density */
            double                      pdf(double r, double p, double q, double x, bool log);                /*!< NegativeBinomial(r,p) probability density */
            double                      lnPdf(double r, double p, double x);                                  /*!< NegativeBinomial(r,p) log_e probability density */
            double                      cdf(double r, double p, double x);                                    /*!< NegativeBinomial(r,p) cumulative probability */
            double                      quantile(double r, double p, double q);                               /*!< NegativeBinomial(r,p) quantile */
            int                         rv(double r, double p, RandomNumberGenerator& rng);                   /*!< NegativeBinomial(r,p) random variable */
	        }
    }
}

#endif
