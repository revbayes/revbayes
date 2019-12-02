/**
 * @file DistributionPoisson
 * This file contains the functions of the poisson distribution.
 *
 * @brief Implementation of the poisson distribution.
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


#ifndef DistributionPoisson_H
#define DistributionPoisson_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace Poisson {
        
            double                      pdf(double lambda, int x);                                            /*!< Poisson(lambda) probability */
            double                      lnPdf(double lambda, int x);                                          /*!< Log of the Poisson(lambda) probability */
            double                      cdf(double lambda, int x);                                            /*!< Poisson(lambda) cumulative probability */
            double                      quantile(double lambda, double p);                                    /*!< Poisson(lambda) quantile */
            int                         rv(double lambda, RandomNumberGenerator& rng);                        /*!< Poisson(lambda) random variable */
        }
    }
}

#endif
