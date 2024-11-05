/**
 * @file DistributionMultinomial
 * This file contains the functions of the multinomial distribution.
 *
 * @brief Implementation of the multinomial distribution.
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


#ifndef DistributionMultinomial_H
#define DistributionMultinomial_H

#include <cstdlib>
#include <vector>

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Multinomial {
            double                      pdf(const std::vector<double>& p, const std::vector<unsigned int>& x);      /*!< Multinomial(p[]) probability */
            double                      pdf(const std::vector<double>& p, const std::vector<std::int64_t>& x);              /*!< Multinomial(p[]) probability */
            double                      pdf(const std::vector<double>& p, const std::vector<double>& x);            /*!< Multinomial(p[]) probability */
            double                      lnPdf(const std::vector<double>& p, const std::vector<unsigned int>& x);    /*!< Multinomial(p[]) log_e probability */
            double                      lnPdf(const std::vector<double>& p, const std::vector<std::int64_t>& x);            /*!< Multinomial(p[]) log_e probability */
            double                      lnPdf(const std::vector<double>& p, const std::vector<double>& x);          /*!< Multinomial(p[]) log_e probability */
            std::vector<std::int64_t>           rv(const std::vector<double>& p, RandomNumberGenerator& rng);               /*!< Multinomial(p[]) random variable */
            std::vector<std::int64_t>           rv(const std::vector<double>& p, size_t n, RandomNumberGenerator& rng);     /*!< Multinomial(p[]) random variable */
        }
    }
}

#endif
