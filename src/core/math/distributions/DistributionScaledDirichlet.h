/**
 * @file DistributionDirichlet
 * This file contains the functions of the dirichlet distribution.
 *
 * @brief Implementation of the dirichlet distribution.
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


#ifndef DistributionDirichlet_H
#define DistributionDirichlet_H

#include <vector>
#include "RbVector.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;


    namespace RbStatistics {
    
        namespace ScaledDirichlet {
            double                      pdf(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& z);        /*!< Dirichlet(a[]) probability density */
            double                      lnPdf(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& z);      /*!< Dirichlet(a[]) log_e probability density */
            RbVector<double>            rv(const std::vector<double>& a, const std::vector<double>& b, RandomNumberGenerator& rng);           /*!< Dirichlet(a[]) random variable */
        }
    }
}

#endif
