/* 
 * File:   DistributionInverseWishart.h
 * Author: JPH, RT
 *
 * Created on February 5, 2015
 */

#ifndef DISTRIBUTIONDECOMPOSEDINVERSEWISHART_H
#define	DISTRIBUTIONDECOMPOSEDINVERSEWISHART_H

#include <cstddef>

#include "MatrixReal.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    namespace RbStatistics {
        
        namespace DecomposedInverseWishart {
            
            // r: matrix of correlation coefficients
            // nu: degrees of freedom
            double                      pdf(double nu, const MatrixReal& r);                  /*!< Decomposed inverse Wishart probability density */
            double                      lnPdf(double nu, const MatrixReal& r);                /*!< Decomposed inverse Wishart log_e probability density */
            MatrixReal                  rv(size_t k, double nu, RandomNumberGenerator& rng);            /*!< Decomposed inverse Wishart random variable */
        }
    }
}


#endif

