/**
 * @file DistributionBinomial
 * This file contains the functions of the binomial distribution.
 *
 * @brief Implementation of the binomial distribution.
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


#ifndef DistributionBinomial_H
#define DistributionBinomial_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Binomial {
        
            double                      pdf(double n, double p, double x);                                    /*!< Binomial(n,p) probability density */
            double                      pdf(double n, double p, double q, double x, bool log);                /*!< Binomial(n,p) probability density */
            double                      lnPdf(double n, double p, double x);                                  /*!< Binomial(n,p) log_e probability density */
            double                      cdf(double n, double p, double x);                                    /*!< Binomial(n,p) cumulative probability */
            double                      quantile(double n, double p, double q);                               /*!< Binomial(n,p) quantile */
            int                         rv(double n, double p, RandomNumberGenerator& rng);             /*!< Binomial(n,p) random variable */
	
            double                      do_search(double y, double *z, double p, double n, double pr, double incr);
        }
    }
}

#endif
