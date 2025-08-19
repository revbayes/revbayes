/**
 * @file
 * This file contains the declaration of RateMatrix_Empirical, which is
 * class that holds a rate matrix for a continuous-time Markov model.
 *
 * @brief Declaration of RateMatrix_Empirical
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-12-05 12:20:37 -0800 (Wed, 05 Dec 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix.h 1901 2012-12-05 20:20:37Z boussau $
 */

#ifndef RateMatrix_Empirical_H
#define RateMatrix_Empirical_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;

    
    class RateMatrix_Empirical : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Empirical(size_t n);                                                                                                                             //!< Construct rate matrix with n states
        RateMatrix_Empirical(const RateMatrix_Empirical& m);                                                                                                        //!< Copy constructor
        virtual                             ~RateMatrix_Empirical(void);                                                                                            //!< Destructor
        
        // overloaded operators
        RateMatrix_Empirical&               operator=(const RateMatrix_Empirical& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_Empirical*               clone(void) const;
        void                                update(void);
        
    protected:
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        
    };
    
}

#endif

