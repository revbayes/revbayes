/**
 * @file
 * This file contains commonly used statistics functions that are used
 * in RevBayes. The probability density (pdf), log of the probability
 * density (lnPdf), cumulative probability (cdf), and quantiles of
 * common probability distributions.
 *
 * @brief Namespace containing statistical functions
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */
#ifndef RbStatisticsHelper_H
#define RbStatisticsHelper_H

#include <cstdint>
#include <vector>

#include "RandomNumberGenerator.h"
#include "RbException.h"

namespace RevBayesCore {

    namespace RbStatistics {
    
        namespace Helper {
        
            double                      dppConcParamFromNumTables(double tables, double num);
            double                      dppExpectNumTableFromConcParam(double conp, double num);
            double                      getACT(const std::vector<double> values, std::int64_t* begin = NULL, std::int64_t* end = NULL); //!< Calculate the auto-correlation time for a time series represented as an arbitrary vector of numeric values
            double                      getESS(const std::vector<double> values, std::int64_t* begin = NULL, std::int64_t* end = NULL); //!< Calculate the effective sample size for a time series represented as an arbitrary vector of numeric values
            double                      getRangeMean(const std::vector<double> values, std::int64_t* begin = NULL, std::int64_t* end = NULL);
            double                      getSEM(const std::vector<double> values, std::int64_t* begin = NULL, std::int64_t* end = NULL); //!< Calculate the standard error of the mean of an arbitrary vector of numeric values
            double                      getVarStat(const std::vector<double> values, std::int64_t* begin = NULL, std::int64_t* end = NULL);
            double                      pointChi2(double prob, double v);
            int                         poissonInver(double lambda, RandomNumberGenerator& rng);
            int                         poissonLow(double lambda, RandomNumberGenerator& rng);
            int                         poissonRatioUniforms(double lambda, RandomNumberGenerator& rng);
            double                      rndGamma(double s, RandomNumberGenerator& rng);
            double                      rndGamma1(double s, RandomNumberGenerator& rng);
            double                      rndGamma2(double s, RandomNumberGenerator& rng);
            double                      rndGamma3(double a, RandomNumberGenerator& rng);
            double                      rndGamma4(double s, RandomNumberGenerator& rng);
		
            template <class T> void		randomlySelectFromVectorWithReplacement(std::vector<T>& sourceV, std::vector<T>& destV, size_t k, RandomNumberGenerator& rng) {
            
                destV.clear();
                for (size_t i=0; i<k; i++)
                    destV.push_back( sourceV[(int)(rng.uniform01()*(sourceV.size()))] );
            }
            template <class T> void		randomlySelectFromVectorWithoutReplacement(std::vector<T>& sourceV, std::vector<T>& destV, size_t k, RandomNumberGenerator& rng) {
            
                if ( (size_t)sourceV.size() < k )
                    throw (RbException("Attempting to sample too many elements from source vector"));
                destV.clear();
                std::vector<T> tmpV = sourceV;
                size_t n = tmpV.size();
                for (size_t i=0; i<k; i++) {
                    size_t whichElement = (size_t)(rng.uniform01()*(n-i));
                    destV.push_back( tmpV[whichElement] );
                    tmpV[whichElement] = tmpV[n-i-1];
                }
            }
            template <class T> void		permuteVector(std::vector<T>& v, RandomNumberGenerator* rng) {
            
                RbStatistics::Helper::randomlySelectFromVectorWithoutReplacement<T>(v, v, v.size(), rng);
            }
        }
        
	}
}


#endif
