//
//  DistributionVarianceGamma.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/8/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef DistributionVarianceGamma_h
#define DistributionVarianceGamma_h

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    namespace RbStatistics {
        
        namespace VarianceGamma {
            
            double                      pdf(double mu, double kappa, double tau, double t, double x);                                            /*!< Poisson(lambda) probability */
            double                      lnPdf(double mu, double kappa, double tau, double t, double x);                                          /*!< Log of the Poisson(lambda) probability */
            double                      cdf(double mu, double kappa, double tau, double t, double x);                                            /*!< Poisson(lambda) cumulative probability */
            double                      quantile(double mu, double kappa, double tau, double t, double p);                                    /*!< Poisson(lambda) quantile */
            double                      rv(double mu, double kappa, double tau, double t, RandomNumberGenerator& rng);                        /*!< Poisson(lambda) random variable */
        }
    }
}

#endif /* DistributionVarianceGamma_h */
