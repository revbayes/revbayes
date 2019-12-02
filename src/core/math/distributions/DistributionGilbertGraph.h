//
//  DistributionGilbertGraph.h
//  revbayes
//
//  Created by Michael R. May on 2017-08-04.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef DistributionGilbertGraph_H
#define DistributionGilbertGraph_H

#include <iostream>
#include "MatrixReal.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    
    namespace RbStatistics {
        
        namespace GilbertGraph {
            
            double                          pdf(double p, const MatrixReal& z);
            double                          lnPdf(double p, const MatrixReal& z);
            MatrixReal                      rv(double p, size_t dim, RandomNumberGenerator& rng);
            
        }
    }
}

#endif /* defined(__revbayes__DistributionGilbertGraph__) */
