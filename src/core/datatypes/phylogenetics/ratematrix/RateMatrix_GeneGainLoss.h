/**
 * @file
 * This file contains the declaration of RateMatrix_GeneGainLoss, which is a
 * class that holds a rate matrix for the chromosome number evolution model.
 *
 * @brief Declaration of RateMatrix_GeneGainLoss
 *
 * (c) copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef RateMatrix_GeneGainLoss_H
#define RateMatrix_GeneGainLoss_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_GeneGainLoss : public AbstractRateMatrix {
        
    public:
        RateMatrix_GeneGainLoss(size_t n, AbstractRateMatrix::METHOD m);                                                  //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_GeneGainLoss(void);                     //!< Destructor
        
        // RateMatrix functions
        double                          averageRate(void) const;
        RateMatrix_GeneGainLoss*        clone(void) const;
        void                            updateInternalRateMatrix(void);
        void                            setBirth(double b);
        void                            setDeath(double d);
        
    private:
        double                          birth;
        double                          death;

        void                            buildRateMatrix(void);
    };
    
}

#endif /* defined(RateMatrix_GeneGainLoss_H) */
