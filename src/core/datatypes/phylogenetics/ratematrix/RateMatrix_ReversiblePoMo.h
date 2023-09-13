/**
 * @file
 * This file contains the declaration of RateMatrix_ReversiblePoMo, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence. Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio, Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_ReversiblePoMo, a reversible matrix combining polymorphisms and substitutions
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2004-07-03 12:20:37 -0800 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix_ReversiblePoMo.h 1901 2014-07-03 06:15 boussau $
 */


#ifndef RateMatrix_ReversiblePoMo_H
#define RateMatrix_ReversiblePoMo_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"
#include "Simplex.h"
#include "RateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
    class Assignable;

    class RateMatrix_ReversiblePoMo : public AbstractRateMatrix {

    public:

        using RateMatrix::getRate;

        RateMatrix_ReversiblePoMo(size_t n);                                                  //!< Construct rate matrix with n states
        //RateMatrix_ReversiblePoMo(size_t n,  const RateGenerator &Qmut);  //!< Construct rate matrix with n states and a matrix of mutation rates
//        RateMatrix_ReversiblePoMo(size_t n,  const RateGenerator &Qmut, const size_t vps);  //!< Construct rate matrix with n states, a matrix of mutation rates, and a virtual population size
        //RateMatrix_ReversiblePoMo(const size_t n,  const RateGenerator &Qmut, const size_t vps ) ;
        RateMatrix_ReversiblePoMo(const size_t n,  const std::vector<double> &rh, const Simplex p, const size_t vps  );

        RateMatrix_ReversiblePoMo(const RateMatrix_ReversiblePoMo& m) ;

        RateMatrix_ReversiblePoMo&                  operator=(const RateMatrix_ReversiblePoMo &r) ;
        virtual                                    ~RateMatrix_ReversiblePoMo(void);                     //!< Destructor

        // RateMatrix functions
        virtual RateMatrix_ReversiblePoMo&          assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_ReversiblePoMo*                  clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, which are the stationary frequencies of the Q_mut matrix

        void                                        update(void);
        void                                        setExchangeabilityRates(const std::vector<double>& r);
        void                                        setStationaryFrequencies(const Simplex& s);


    private:

        void                                        decomposeState(int state, int &i, int &nt1, int &nt2) ;
        double                                      mutCoeff(int nt1, int nt2) ;
        double                                      computeProbBoundaryMutation(int state1, int state2) ;
        void                                        buildRateMatrix(void) ;
        void                                        calculateCijk(void) ;
        void                                        tiProbsEigens(double t, TransitionProbabilityMatrix& P) const ;
        void                                        tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const ;
        void                                        updateEigenSystem(void);

        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const ;
        inline void                                 squareMatrix( TransitionProbabilityMatrix& P, TransitionProbabilityMatrix& P2) const ;
        
        
        std::vector<double>                         rho;                                        //!< Holds the exchangeabilities
        Simplex                                     pi;                                         //!< Holds the stationary frequencies
        
        size_t                                      N;                                            //!< Number of individuals in idealized population
        size_t                                      matrixSize;                                 //!< Number of elements in a row or column of the rate matrix
        double                                      precision;                                  //!< Precision for exponentiation through repeated squaring
        
        
    };

}

#endif /* defined(__RateMatrix_ReversiblePoMo__) */
