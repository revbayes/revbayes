//
//  DistributionT.hpp
//  revbayes_work
//
//  Created by Dismukes, Wade T [EEOBS] on 10/26/16.
//  Copyright © 2016 Dismukes, Wade T [EEOBS]. All rights reserved.
//

/**
 * @file DistributionF
 * This file contains the functions of the F distribution.
 *
 * @brief Implementation of the F distribution.
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


#ifndef DistributionT_H
#define DistributionT_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    namespace RbStatistics {
        
        namespace StudentT {
            
            double                      pdf(double v, double x);                                    /*!< T(v) probability density */
            double                      pdf(double v, double x, bool give_log);                     /*!< T(v) probability density */
            double                      lnPdf(double v, double x);                                  /*!< T(v) log_e probability density */
            double                      cdf(double v, double x);                                /*!< T(v) cumulative probability */
            double                      quantile(double v, double p);                               /*!< T(v) quantile */
            double                      rv(double v, RandomNumberGenerator& rng);                   /*!< T(v) random variable */
            
        }
    }
}

#endif
