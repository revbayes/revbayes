/**
 * @file
 * This file contains the declaration of RateMatrix_PoMoBalanceKN, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence and balancing selection. This  is an extension
 * of the model Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio,
 * Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic
 * models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_PoMoBalanceKN, a non-reversible matrix combining polymorphisms and substitutions with balancing selection
 *
 * (c) Copyright 2023-
 * @date Last modified: $Date: 2023-12-12 17:32:37 -0800 $
 * @author The RevBayes Development Core Team (Svitlana Braichenko)
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMoBalanceKN.h 1901 2022-08-08 06:15 boussau $
 */


#ifndef RateMatrix_PoMoBalanceKN_H
#define RateMatrix_PoMoBalanceKN_H

#include <stddef.h>
#include <vector>

#include "AbstractRateMatrix.h"
#include "Simplex.h"
#include "RateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
    class Assignable;

    class RateMatrix_PoMoBalanceKN : public AbstractRateMatrix {

    public:

        using RateMatrix::getRate;

        //RateMatrix_PoMoBalanceKN(size_t num_states) ;
        RateMatrix_PoMoBalanceKN( long num_states, long in_k, long in_n, long in_nmr )  ;
        RateMatrix_PoMoBalanceKN(const RateMatrix_PoMoBalanceKN& m) ;

        RateMatrix_PoMoBalanceKN&                         operator=(const RateMatrix_PoMoBalanceKN &r) ;
        virtual                                    ~RateMatrix_PoMoBalanceKN(void);                     //!< Destructor

        // RateMatrix functions
        virtual RateMatrix_PoMoBalanceKN&           assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMoBalanceKN*                   clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, which are the stationary frequencies of the Q_mut matrix

        void                                        update(void);
        void                                        setK( long &na );
        void                                        setN( long &ni );
        void                                        setMu(  const std::vector<double> &m );
        void                                        setPhi( const std::vector<double> &f );
        void                                        setBeta( const std::vector<double> &b );
        void                                        setB( const std::vector<long> &Bf );


    private:
        void                                        buildRateMatrix(void) ;
        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const ;

        long                                        K;
        long                                        N;
        std::vector<double>                         mu;   
        std::vector<double>                         phi;
        std::vector<double>                         beta;                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<long>                           B;
        // std::vector<double>                         stationaryVector;                    //!< Holds the stationary frequencies

    };

}

#endif /* defined(__RateMatrix_PoMoBalanceKN__) */
