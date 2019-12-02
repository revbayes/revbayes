//
//  DistributionCompoundPoissonNormal.h
//  revbayes-proj
//
//  Created by Michael Landis on 11/19/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef DistributionCompoundPoissonNormal_H
#define DistributionCompoundPoissonNormal_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    namespace RbStatistics {
        
        namespace CompoundPoissonNormal {
            
            double                      pdf(double lambda, double mu, double sigma, double x);                                            /*!< Poisson(lambda) probability */
            double                      lnPdf(double lambda, double mu, double sigma, double x);                                          /*!< Log of the Poisson(lambda) probability */
            double                      cdf(double lambda, double mu, double sigma, double x);                                            /*!< Poisson(lambda) cumulative probability */
            double                      quantile(double lambda, double mu, double sigma, double p);                                    /*!< Poisson(lambda) quantile */
            double                      rv(double lambda, double mu, double sigma, RandomNumberGenerator& rng);                        /*!< Poisson(lambda) random variable */
        }
    }
}

#endif /* defined(__revbayes_proj__DistributionCompoundPoissonNormal__) */
