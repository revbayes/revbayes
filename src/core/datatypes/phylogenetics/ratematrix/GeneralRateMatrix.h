#ifndef GeneralRateMatrix_H
#define GeneralRateMatrix_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
class Assignable;
    
    
    /**
     * \brief Abstract class for general rate matrices.
     *
     * This is an abstract class for general rate matrices that implementes the RateMatrix interface.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2013-04-25
     *
     */
    class GeneralRateMatrix : public AbstractRateMatrix {
        
    public:
        GeneralRateMatrix(size_t n, bool rescale_to_one=true);                                                                            //!< Construct rate matrix with n states
        virtual                            ~GeneralRateMatrix(void);                                                                    //!< Destructor
        
        
        // public methods
        
        // pure virtual methods you have to overwrite
        virtual GeneralRateMatrix&          assign(const Assignable &m);
        virtual double                      averageRate(void) const;                                                                    //!< Calculate the average rate
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const = 0;   //!< Calculate the transition matrix
        virtual GeneralRateMatrix*          clone(void) const = 0;
        const std::vector<double>&          getTransitionRates(void) const;
        virtual std::vector<double>         getStationaryFrequencies(void) const;                                                       //!< Return the stationary frequencies
        bool                                isTimeReversible(void);                                                                     //!< Return whether the rate matrix is time reversible
        void                                setTransitionRates(const std::vector<double> &tr);
//        void                                setStationaryFrequencies(const std::vector<double>& f);                                     //!< Directly set the stationary frequencies
        virtual void                        update(void);                                                                               //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)

    protected:

        // members
//        std::vector<double>                 stationary_freqs;
        std::vector<double>                 transition_rates;
        bool                                rescale_to_one;
        
    };
    
    
}

#endif


