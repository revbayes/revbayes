/**
 * @file
 * This file contains the declaration of RateMatrix_ChromosomesPloidy, which is a
 * class that holds a rate matrix for the chromosome number evolution model. 
 *
 * @brief Declaration of RateMatrix_ChromosomesPloidy
 *
 * (c) copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef __RateMatrix_ChromosomesPloidy__
#define __RateMatrix_ChromosomesPloidy__

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_ChromosomesPloidy : public AbstractRateMatrix {
        
    public:
        RateMatrix_ChromosomesPloidy(size_t n);                                                  //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_ChromosomesPloidy(void);                     //!< Destructor
        
        // RateMatrix functions
        double                          averageRate(void) const;
        void                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_ChromosomesPloidy*   clone(void) const;
        std::vector<double>             getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in this model I don't know them
        void                            update(void);
        void                            setGamma_d(double g);
        void                            setDelta_d(double d);
        void                            setRho_d(double r);
        void                            setEta_d(double e);
        void                            setGamma_dl(double l);
        void                            setDelta_dl(double d);
        void                            setGamma_p(double g);
        void                            setDelta_p(double d);
        void                            setRho_p(double r);
        void                            setEta_p(double e);
        void                            setGamma_pl(double l);
        void                            setDelta_pl(double d);
        void                            setBeta(double b);
        
        
    private:
        double                          gamma_d;
        double                          delta_d;
        double                          rho_d;
        double                          eta_d;
        double                          gamma_dl;
        double                          delta_dl;
        double                          gamma_p;
        double                          delta_p;
        double                          rho_p;
        double                          eta_p;
        double                          gamma_pl;
        double                          delta_pl;
        double                          beta;
        size_t                          matrix_size;                         //!< Number of elements in a row or column of the rate matrix
        size_t                          max_chromo;                         
        std::vector<double>             stationary_freqs;                    //!< Holds the stationary frequencies

        void                            buildRateMatrix(void);
    };
    
}

#endif /* defined(__RateMatrix_ChromosomesPloidy__) */
