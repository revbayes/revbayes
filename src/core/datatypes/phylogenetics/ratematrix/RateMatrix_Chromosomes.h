/**
 * @file
 * This file contains the declaration of RateMatrix_Chromosomes, which is a
 * class that holds a rate matrix for the chromosome number evolution model. 
 *
 * @brief Declaration of RateMatrix_Chromosomes
 *
 * (c) copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef __RateMatrix_Chromosomes__
#define __RateMatrix_Chromosomes__

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_Chromosomes : public AbstractRateMatrix {
        
    public:
        RateMatrix_Chromosomes(size_t n);                                                  //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_Chromosomes(void);                     //!< Destructor
        
        // RateMatrix functions
        double                          averageRate(void) const;
        void                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_Chromosomes*         clone(void) const;
        std::vector<double>             getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in this model I don't know them
        void                            update(void);
        void                            setGamma(double g);
        void                            setDelta(double d);
        void                            setRho(double r);
        void                            setEta(double e);
        void                            setGamma_l(double l);
        void                            setDelta_l(double d);
        
        
    private:
        double                          gamma;
        double                          delta;
        double                          rho;
        double                          eta;
        double                          gamma_l;
        double                          delta_l;
        size_t                          matrixSize;                         //!< Number of elements in a row or column of the rate matrix
        std::vector<double>             stationary_freqs;                    //!< Holds the stationary frequencies

        void                            buildRateMatrix(void);
    };
    
}

#endif /* defined(__RateMatrix_Chromosomes__) */
