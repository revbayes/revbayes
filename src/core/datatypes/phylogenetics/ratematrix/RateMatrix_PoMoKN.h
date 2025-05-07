/**
 * @file
 * This file contains the declaration of RateMatrix_PoMoKN, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence. Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio, Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_PoMoKN, a reversible matrix combining polymorphisms and substitutions
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2004-07-03 12:20:37 -0800 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMoKN.h 1901 2014-07-03 06:15 boussau $
 */


#ifndef RateMatrix_PoMoKN_H
#define RateMatrix_PoMoKN_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"
#include "Simplex.h"
#include "RateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
    class Assignable;

    class RateMatrix_PoMoKN : public AbstractRateMatrix {

    public:

        using RateMatrix::getRate;

        //RateMatrix_PoMoKN(size_t num_states) ;
        RateMatrix_PoMoKN(long num_states, long in_k, long in_n, long in_nmr)  ;
        RateMatrix_PoMoKN(const RateMatrix_PoMoKN& m) ;

        RateMatrix_PoMoKN&                         operator=(const RateMatrix_PoMoKN &r) ;
        virtual                                    ~RateMatrix_PoMoKN(void);                     //!< Destructor

        // RateMatrix functions
        virtual RateMatrix_PoMoKN&                  assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMoKN*                          clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, which are the stationary frequencies of the Q_mut matrix

        void                                        update(void);
        void                                        setK( long &na );
        void                                        setN( long &ni );
        void                                        setMu(  const std::vector<double> &m );
        void                                        setPhi( const std::vector<double> &f );


    private:
        void                                        buildRateMatrix(void) ;
        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const ;
        
        long                                        K;
        long                                        N;
        std::vector<double>                         mu;   
        std::vector<double>                         phi;    
        std::vector<double>                         stationaryVector;                    //!< Holds the stationary frequencies

    };

}

#endif /* defined(__RateMatrix_PoMoKN__) */
