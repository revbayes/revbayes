/**
 * @file DistributionGeometric
 * This file contains the functions of the geometric distribution.
 *
 * @brief Implementation of the geometric distribution.
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


#ifndef DistributionGeometric_H
#define DistributionGeometric_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace Geometric {
        
            double                      pdf(int n, double p);                                   /*!< Geometric(p) probability density */
            double                      pdf(int n, double p, bool asLog);                       /*!< Geometric(p) probability density */
            double                      lnPdf(int n, double p);                                 /*!< Geometric(p) log_e probability density */
            double                      cdf(int n, double p);                                   /*!< Geometric(p) cumulative probability */
            int                         quantile(double q, double p);                           /*!< Geometric(p) quantile */
            int                         rv(double p, RandomNumberGenerator& rng);               /*!< Geometric(p) random variable */
        }
    }
}

#endif
