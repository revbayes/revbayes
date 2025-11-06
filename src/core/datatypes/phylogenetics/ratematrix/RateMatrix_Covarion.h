//
//  RateMatrix_Covarion.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/4/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef RateMatrix_Covarion_hpp
#define RateMatrix_Covarion_hpp

#include "GeneralRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    /**
     * @brief Free rate matrix class.
     *
     * This class implements the free rate matrix .
     * The resulting rate matrix is computed by:
     *
     *      |   -            r[1]      r[2]    ...    r[k]    |
     *      |                                                 |
     *      |  r[k+1]         -       r[k+2]   ...   r[2k]    |
     * Q =  |                                                 |
     *      |  r[(k-2)k+1]          ...      -     r[(k-1)k]  |
     *      |                                                 |
     *      |  r[(k-1)k+1]          ...                -      |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-07-04, version 1.0
     */
    class RateMatrix_Covarion : public GeneralRateMatrix {
        
    public:
        RateMatrix_Covarion(size_t n, size_t k);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_Covarion(size_t n, size_t k, bool r);
        RateMatrix_Covarion(const RateMatrix_Covarion& m);                                                                                //!< Copy constructor
        virtual                         ~RateMatrix_Covarion(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_Covarion&                operator=(const RateMatrix_Covarion& r);
        std::vector<int>                    get_emitted_letters(void) const;
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_Covarion*                clone(void) const;
        void                                fillRateMatrix(void);
        void                                setRateMatrices(const RbVector<MatrixReal> &rm);
        void                                setSwitchRates(const MatrixReal &sr);
        void                                setClockRates(const RbVector<double>& cr);                                     //!< Directly set the stationary
        void                                update(void);
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        bool                                rescale;
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        // matrix member variables
        RbVector<MatrixReal>                rate_matrices;
        MatrixReal                          switch_rates;
        RbVector<double>                    clock_rates;
        size_t                              num_classes;
        size_t                              num_states_per_class;
        
        
    };
    
}


#endif /* RateMatrix_Covarion_hpp */
