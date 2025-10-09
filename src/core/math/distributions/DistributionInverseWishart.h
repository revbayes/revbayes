/* 
 * File:   DistributionInverseWishart.h
 * Author: nl
 *
 * Created on 15 juillet 2014, 10:13
 */

#ifndef DISTRIBUTIONINVERSEWISHART_H
#define	DISTRIBUTIONINVERSEWISHART_H


#include <cstddef>
#include <vector>

#include "MatrixReal.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    
    namespace RbStatistics {
        
        namespace InverseWishart {
            
            // standard form:
            // sigma0: scale matrix
            // df: degrees of freedom
            // by definition: Z = \sum_i X_i X_i', where X_i ~ Normal(0, covariance = sigma0)
            double                      pdf(const MatrixReal& sigma0, size_t df, const MatrixReal& z);                  /*!< Dirichlet(a[]) probability density */
            double                      lnPdf(const MatrixReal& sigma0, size_t df, const MatrixReal& z);                /*!< Dirichlet(a[]) log_e probability density */
            MatrixReal                  rv(const MatrixReal& sigma0, size_t df, RandomNumberGenerator& rng);            /*!< Dirichlet(a[]) random variable */
            MatrixReal                  rvCovariance(const MatrixReal& sigma0, size_t df, RandomNumberGenerator& rng);            /*!< Dirichlet(a[]) random variable */

            // InverseWishart of parameters sigma0 = Diagonal(kappa) and df degrees of freedom
            double                      pdf(const std::vector<double>& kappa, size_t df, const MatrixReal& z);          /*!< Dirichlet(a[]) probability density */
            double                      lnPdf(const std::vector<double>& kappa, size_t df, const MatrixReal& z);        /*!< Dirichlet(a[]) log_e probability density */
            MatrixReal                  rv(const std::vector<double>& kappa, size_t df, RandomNumberGenerator& rng);    /*!< Dirichlet(a[]) random variable */

            // InverseWishart of parameters sigma0 = kappa * I and df degrees of freedom
            double                      pdf(double kappa, size_t df, const MatrixReal& z);                              /*!< Dirichlet(a[]) probability density */
            double                      lnPdf(double kappa, size_t df, const MatrixReal& z);                            /*!< Dirichlet(a[]) log_e probability density */
            MatrixReal                  rv(double kappa, size_t dim, size_t df, RandomNumberGenerator& rng);            /*!< Dirichlet(a[]) random variable */

        }
    }
}


#endif	/* DISTRIBUTIONINVERSEWISHART_H */

