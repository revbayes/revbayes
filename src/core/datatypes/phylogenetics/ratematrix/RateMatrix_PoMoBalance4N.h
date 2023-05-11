/**
 * @file
 * This file contains the declaration of RateMatrix_PoMoBalance4N, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence and balancing selection. This  is an extension
 * of the model Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio,
 * Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic
 * models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_PoMoBalance4N, a reversible matrix combining polymorphisms and substitutions with balancing selection
 *
 * (c) Copyright 2022-
 * @date Last modified: $Date: 2022-08-08 12:20:37 -0800 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMoBalance4N.h 1901 2022-08-08 06:15 boussau $
 */


#ifndef RateMatrix_PoMoBalance4N_H
#define RateMatrix_PoMoBalance4N_H

#include <stddef.h>
#include <vector>

#include "AbstractRateMatrix.h"
#include "Simplex.h"
#include "RateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
    class Assignable;

    class RateMatrix_PoMoBalance4N : public AbstractRateMatrix {

    public:

        using RateMatrix::getRate;

        //RateMatrix_PoMoBalance4N(size_t num_states) ;
        RateMatrix_PoMoBalance4N( long num_states, long in_n )  ;
        RateMatrix_PoMoBalance4N(const RateMatrix_PoMoBalance4N& m) ;

        RateMatrix_PoMoBalance4N&                         operator=(const RateMatrix_PoMoBalance4N &r) ;
        virtual                                    ~RateMatrix_PoMoBalance4N(void);                     //!< Destructor

        // RateMatrix functions
        virtual RateMatrix_PoMoBalance4N&           assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMoBalance4N*                   clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, which are the stationary frequencies of the Q_mut matrix

        void                                        update(void);
        void                                        setN( long &ni );
        void                                        setMu(  const std::vector<double> &m );
        void                                        setPhi( const std::vector<double> &f );
        void                                        setBeta( const std::vector<double> &b );
        void                                        setB( const std::vector<long> &Bf );


    private:
        void                                        buildRateMatrix(void) ;
        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const ;
        
        long                                        N;
        std::vector<double>                         mu;   
        std::vector<double>                         phi;
        std::vector<double>                         beta;                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<long>                           B;
        // std::vector<double>                         stationaryVector;                    //!< Holds the stationary frequencies

    };

}

#endif /* defined(__RateMatrix_PoMoBalance4N__) */
