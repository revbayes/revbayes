/**
 * @file
 * This file contains the declaration of RateMatrix_PoMo4, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence. Presented in De Maio, Schlötterer and Kosiol, Molecular Biology and Evolution 30:10 2249-62 (2013).
 * Parameters:
 * N: Number of individuals in the idealized populations.
 * mu: matrix of ordered mutation rates of size 4*4: 12 mutation rates: AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG  describing mutations between monomorphic states
 * In the original publication, the model is symmetric: unordered mutation rates of size 6: AC, AG, AT, CG, CT, GT, CA=AC, GA=AG, etc...
 * s: vector of fitness parameters of size 4: sA, sC, sG, sT
 // * pi: vector of root frequencies of size 4: piA, piC, piG, piT
 *
 * @brief Declaration of RateMatrix_PoMo4, a matrix combining polymorphisms and substitutions
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2004-07-03 12:20:37 -0800 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMo4.h 1901 2014-07-03 06:15 boussau $
 */


#ifndef RateMatrix_PoMo4_H
#define RateMatrix_PoMo4_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"
#include "RateGenerator.h"
#include "RateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_PoMo4 : public AbstractRateMatrix {
        
    public:
        
        using RateMatrix::getRate;
        
        RateMatrix_PoMo4(size_t n);                                                  //!< Construct rate matrix with n states
        RateMatrix_PoMo4(size_t n, const size_t vps, const std::vector<double> &mr, const std::vector<double> &sc);  //!< Construct rate matrix with n states, a vector of mutation rates, and a vector of selection coefficients
        RateMatrix_PoMo4(size_t n, const size_t vps, const RateGenerator &mm, const std::vector<double> sc);  //!< Construct rate matrix with n states, a matrix of mutation rates, and a vector of selection coefficients
        
        virtual                         ~RateMatrix_PoMo4(void);                     //!< Destructor
        
        // RateMatrix functions
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMo4*                           clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in the PoMo4 model I don't know them
        
        void                                        update(void);
        void                                        setMutationRates(const std::vector<double>& mr);
        void                                        setMutationRates(const RateGenerator& mm);
        void                                        setMutationRates(const std::vector<double>& r, const Simplex& s);
        void                                        setSelectionCoefficients(const std::vector<double>& sc);
        
        
    private:
        void                                        buildRateMatrix(void);
        double                                      computeEntryFromMoranProcessWithSelection(size_t state1, size_t state2, double& count1);
        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const;
        inline void                                 squareMatrix( TransitionProbabilityMatrix& P, TransitionProbabilityMatrix& P2) const;
        
        
        size_t                                      N;                          //!< Number of individuals in idealized population
        size_t                                      matrixSize;                 //!< Number of elements in a row or column of the rate matrix
        std::vector < std::vector < double > >      mu;                         //!< Matrix of 12 mutation rates and 0s elsewhere
        std::vector < double >                      s;                          //!< Vector of 4 selection coefficients
        double                                      precision;                  //!< Precision for exponentiation through repeated squaring
        std::vector<double>                         stationary_freqs;           //!< Holds the stationary frequencies
        
    };
    
}

#endif /* defined(__RateMatrix_PoMo4__) */
