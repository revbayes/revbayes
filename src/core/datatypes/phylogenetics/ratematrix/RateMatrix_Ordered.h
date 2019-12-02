/**
 * @file
 * This file contains the declaration of RateMatrix_Ordered, which is a
 * class that holds a rate matrix for the chromosome number evolution model.
 *
 * @brief Declaration of RateMatrix_Ordered
 *
 * (c) copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef RateMatrix_Ordered_H
#define RateMatrix_Ordered_H

#include "AbstractRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_Ordered : public AbstractRateMatrix {
        
    public:
        RateMatrix_Ordered(size_t n);                                                  //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_Ordered(void);                     //!< Destructor
        
        // RateMatrix functions
        double                          averageRate(void) const;
        void                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_Ordered*             clone(void) const;
        std::vector<double>             getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in this model I don't know them
        void                            update(void);
        void                            setLambda(double l);
        void                            setMu(double m);
        
        
    private:
        double                          lambda;
        double                          mu;
        size_t                          matrix_size;                                    //!< Number of elements in a row or column of the rate matrix
        std::vector<double>             stationary_freqs;                               //!< Holds the stationary frequencies
        
        void                            buildRateMatrix(void);
        void                            exponentiateMatrixByScalingAndSquaring(double t,  TransitionProbabilityMatrix& p) const;
        inline void                     multiplyMatrices(TransitionProbabilityMatrix& p,  TransitionProbabilityMatrix& q,  TransitionProbabilityMatrix& r) const;

    };
    
}

#endif
