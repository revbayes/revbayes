/**
 * @file DistributionBivariatePoisson
 * This file contains the functions of the bivariate poisson distribution.
 *
 * @brief Implementation of the bivariate poisson distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Alexander E. Zarebski
 * @license GPL version 3
 * @version 1.0
 * @since 2018-12-19, version 1.0
 *
 * $Id$
 */


#ifndef DistributionBivariatePoisson_H
#define DistributionBivariatePoisson_H

#include <vector>

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace BivariatePoisson {
        
          double                      lnPdf(double theta1, double theta2, double theta0, std::vector<std::int64_t> x);
          std::vector<std::int64_t>            rv(double theta1, double theta2, double theta0, RandomNumberGenerator& rng);
        }
    }
}

#endif
