//
//  DistributionWishart.h
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-28.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef DistributionWishart_H
#define DistributionWishart_H

#include <iostream>


#include "MatrixReal.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    
    namespace RbStatistics {
        
        namespace Wishart {
            
            // standard form:
            // omega0: scale matrix
            // df: degrees of freedom
            // by definition: Z = \sum_i X_i X_i', where X_i ~ Normal(0, covariance = omega0)
            double                          pdf(const MatrixReal& omega0, size_t df, const MatrixReal& z);                  /*!< Dirichlet(a[]) probability density */
            double                          lnPdf(const MatrixReal& omega0, size_t df, const MatrixReal& z);                /*!< Dirichlet(a[]) log_e probability density */
            MatrixReal                      rv(const MatrixReal& omega0, size_t df, RandomNumberGenerator& rng);            /*!< Dirichlet(a[]) random variable */

            // Wishart of parameters omega0 = kappa * I and df degrees of freedom
            double                          pdf(double kappa, size_t df, const MatrixReal& z);                              /*!< Dirichlet(a[]) probability density */
            double                          lnPdf(double kappa, size_t df, const MatrixReal& z);                            /*!< Dirichlet(a[]) log_e probability density */
            MatrixReal                      rv(double kappa, size_t dim, size_t df, RandomNumberGenerator& rng);            /*!< Dirichlet(a[]) random variable */

        }
    }
}


#endif /* defined(__revbayes__DistributionWishart__) */
