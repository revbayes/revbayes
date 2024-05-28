/**
 * @file
 * This file contains the declaration of RateMatrix_PoMo2N, which is a
 * class that holds a rate matrix combining population-level polymorphisms
 * with inter-species divergence. Presented in Dominik Schrempf, Bui Quang Minh, Nicola De Maio, Arndt von Haeseler, Carolin Kosiol. Reversible polymorphism-aware phylogenetic models and their application to tree inference. Journal of Theoretical Biology, 2016.
 * Parameters:
 * Q_mut: A GTR instantaneous mutation matrix.
 * N: effective population size.
 *
 * @brief Declaration of RateMatrix_PoMo2N, a reversible matrix combining polymorphisms and substitutions
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2004-07-03 12:20:37 -0800 $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix_PoMo2N.h 1901 2014-07-03 06:15 boussau $
 */


#ifndef RateMatrix_PoMo2N_H
#define RateMatrix_PoMo2N_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"
#include "Simplex.h"
#include "RateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
    class Assignable;

    class RateMatrix_PoMo2N : public AbstractRateMatrix {

    public:

        using RateMatrix::getRate;

        //RateMatrix_PoMo2N(size_t num_states) ;
        RateMatrix_PoMo2N(long num_states, long in_n )  ;
        RateMatrix_PoMo2N(const RateMatrix_PoMo2N& m) ;

        RateMatrix_PoMo2N&                         operator=(const RateMatrix_PoMo2N &r) ;
        virtual                                    ~RateMatrix_PoMo2N(void);                     //!< Destructor

        // RateMatrix functions
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMo2N*                          clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, which are the stationary frequencies of the Q_mut matrix

        void                                        update(void);
        void                                        setN( long &ni );
        void                                        setMu(  const std::vector<double> &m );
        void                                        setPhi( const std::vector<double> &f );


    private:
        void                                        buildRateMatrix(void) ;
        void                                        computeExponentialMatrixByRepeatedSquaring(double t, TransitionProbabilityMatrix& P ) const ;
        
        long                                        N;
        std::vector<double>                         mu;   
        std::vector<double>                         phi;    
        std::vector<double>                         stationaryVector;                    //!< Holds the stationary frequencies

    };

}

#endif /* defined(__RateMatrix_PoMo2N__) */
