//
//  DistributionInverseInverseGamma.h
//  rb_mlandis
//
//  Created by Michael Landis on 3/19/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__DistributionInverseInverseGamma__
#define __rb_mlandis__DistributionInverseInverseGamma__

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {

        namespace InverseGamma {

            double                      pdf(double a, double b, double x);                                    /*!< InverseGamma(a,b) probability density */
            double                      pdf(double a, double b, double x, bool isLog);                        /*!< InverseGamma(a,b) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< Log of the InverseGamma(a,b) probability density */
            double                      cdf(double a, double b, double x);                                    /*!< InverseGamma(a,b) cumulative probability */
            double                      quantile(double a, double b, double p);                               /*!< InverseGamma(a,b) quantile */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< InverseGamma(a,b) random variable */

        }
    }
}

#endif /* defined(__rb_mlandis__DistributionInverseInverseGamma__) */
