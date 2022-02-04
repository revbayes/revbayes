#ifndef DistributionLKJ_H
#define DistributionLKJ_H

#include <iostream>
#include "MatrixReal.h"

namespace RevBayesCore {
    
    class RandomNumberGenerator;
    
    
    namespace RbStatistics {
        
        namespace LKJ {
            
            double                          pdf(double eta, const MatrixReal& z);
            double                          lnPdf(double eta, const MatrixReal& z);
            MatrixReal                      rv(double eta, size_t dim, RandomNumberGenerator& rng);

            double                          pdfPartial(double eta, const MatrixReal& z);
            double                          lnPdfPartial(double eta, const MatrixReal& z);
            MatrixReal                      rvPartial(double eta, size_t dim, RandomNumberGenerator& rng);

        }
    }
}

#endif /* defined(DistributionLKJ_H) */
