/**
 * @file
 * This file contains the declaration of RateMatrix_revPoMoBalanceKN, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence and balancing selection. This  is an extension
 * of the model Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio,
 * Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic
 * models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_revPoMoBalanceKN, a reversible matrix combining polymorphisms and substitutions with balancing selection
 *
 * (c) Copyright 2023-
 * @date Last modified: $Date: 2023-12-12 17:32:37 -0800 $
 * @author The RevBayes Development Core Team (Svitlana Braichenko)
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMoBalanceKN.h 1901 2022-08-08 06:15 boussau $
 */
#ifndef RateMatrix_revPoMoBalanceKN_H
#define RateMatrix_revPoMoBalanceKN_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {

    class EigenSystem;
    class TransitionProbabilityMatrix;

    /*
     */

    class RateMatrix_revPoMoBalanceKN : public TimeReversibleRateMatrix {

    public:
        RateMatrix_revPoMoBalanceKN( std::int64_t num_states, std::int64_t in_k, std::int64_t in_n, std::int64_t in_nex );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_revPoMoBalanceKN(const RateMatrix_revPoMoBalanceKN& m);                                                                  //!< Copy constructor
        virtual                             ~RateMatrix_revPoMoBalanceKN(void);                                                              //!< Destructor

        // overloaded operators
        RateMatrix_revPoMoBalanceKN&                            operator=(const RateMatrix_revPoMoBalanceKN& r);

        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMoBalanceKN*                            clone(void) const;

        void                                                    setK( std::int64_t &na );
        void                                                    setN( std::int64_t &ni );
        void                                                    setPi(const std::vector<double> &p );
        void                                                    setRho( const std::vector<double> &r );
        void                                                    setPhi( const std::vector<double> &s );
        void                                                    setBeta( const std::vector<double> &b );
        std::vector<double>                                     getStationaryFrequencies( void ) const;


        void                                                    update(void);

    private:
        void                                                    calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                                    computeOffDiagonal( void );
        void                                                    tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                                    tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                                    updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors

        EigenSystem*                                            eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                                     c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >                      cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case

        std::int64_t                                                    K;
        std::int64_t                                                    N;
        std::vector<double>                                     pi;
        std::vector<double>                                     rho;
        std::vector<double>                                     phi;
        std::vector<double>                                     beta;                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::int64_t>                                        B;
    };

}

#endif
