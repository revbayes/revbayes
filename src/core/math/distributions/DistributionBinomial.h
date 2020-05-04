/**
 * @file DistributionBinomial.h
 * @brief Implementation of the binomial distribution.
 *
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
