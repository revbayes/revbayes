/*
 * File:   DistributionExponentialError.h
 * Author: David Cerny
 *
 * Created on October 8, 2019
 */

#ifndef DISTRIBUTIONEXPONENTIALERROR_H
#define DISTRIBUTIONEXPONENTIALERROR_H


#include <cstddef>
#include <vector>
#include <iostream>

#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    
    namespace RbStatistics {
        
        namespace ExponentialError {
        
            // exponential error of parameters avgDistMat and lambda
            double                      pdf(const AverageDistanceMatrix& avgDistMat, double lambda, const AverageDistanceMatrix& z);   /*!< Exponential error probability density */
            double                      lnPdf(const AverageDistanceMatrix& avgDistMat, double lambda, const AverageDistanceMatrix& z); /*!< Exponential error log_e probability density */
            AverageDistanceMatrix       rv(const AverageDistanceMatrix& avgDistMat, double lambda, RandomNumberGenerator& rng);        /*!< Exponential error random variable */

            // exponential error of parameters distMat = avgDistMat.getDistanceMatrix() and lambda
            double                      pdf(const DistanceMatrix& distMat, double lambda, const AverageDistanceMatrix& z);                    /*!< Exponential error probability density */
            double                      lnPdf(const DistanceMatrix& distMat, double lambda, const AverageDistanceMatrix& z);                  /*!< Exponential error log_e probability density */
            AverageDistanceMatrix       rv(const DistanceMatrix& distMat, double lambda, RandomNumberGenerator& rng);                  /*!< Exponential error random variable */
        
        }
    }
}


#endif    /* DISTRIBUTIONEXPONENTIALERROR_H */
