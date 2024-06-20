#ifndef AbstractRateMatrix_H
#define AbstractRateMatrix_H

#include <cstddef>
#include <vector>

#include "MatrixReal.h"
#include "RateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    /**
     * @brief Abstract rate matrix class.
     *
     * The abstract rate matrix class provides some basic functionality of most
     * rate matrices. It implements some functionality of the interface RateMatrix.
     * The key element of this abstract class is that it hold a matrix of doubles as the internal value
     * which represent the rates of this matrix. Specific behaviour how the rates are updated,
     * how the stationary frequencies are computed, and how the transition probabilities are computed are
     * left for derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-04, version 1.0
     */
    class AbstractRateMatrix : public RateMatrix {
        
    public:
        virtual                            ~AbstractRateMatrix(void);                                                                   //!< Destructor
        
        // overloaded operators
//        std::vector<double>&                operator[](size_t i);                                                                       //!< Subscript operator
//        const std::vector<double>&          operator[](size_t i) const;                                                                 //!< Subscript operator (const)
        
//        std::vector<std::vector<double> >::const_iterator       begin(void) const;
//        std::vector<std::vector<double> >::iterator             begin(void);
//        std::vector<std::vector<double> >::const_iterator       end(void) const;
//        std::vector<std::vector<double> >::iterator             end(void);
        
        // public methods
        double                              getRate(size_t from, size_t to, double rate=1.0) const;
        double                              getRate(size_t from, size_t to, double age, double rate) const;         //!< Calculate the rate from state i to state j over the given time interval scaled by a rate
        void                                rescaleToAverageRate(double r);                                                             //!< Rescale the rate matrix such that the average rate is "r"
        void                                setDiagonal(void);                                                                          //!< Set the diagonal such that each row sums to zero
        virtual std::vector<int>            get_emitted_letters() const;                                                                //!<Find out what alphet letter each state emits        

        // pure virtual methods you have to overwrite
        virtual double                      averageRate(void) const = 0;                                                                //!< Calculate the average rate
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const = 0;   //!< Calculate the transition matrix
        virtual void                        calculateTransitionProbabilitiesForStochasticMapping(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        virtual AbstractRateMatrix*         clone(void) const = 0;
        virtual std::vector<double>         getStationaryFrequencies(void) const = 0;                                                   //!< Return the stationary frequencies
        MatrixReal                          getRateMatrix(void) const;
        virtual void                        update(void) = 0;                                                                           //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)
        virtual MatrixReal                  getStochasticMatrix(size_t n) const;
        virtual double                      getDominatingRate(void) const;
        virtual bool                        simulateStochasticMapping(double startAge, double endAge, double rate,std::vector<size_t>& transition_states, std::vector<double>& transition_times) const;
        

    protected:
        // prevent instantiation
        AbstractRateMatrix(size_t n);                                                                                                   //!< Construct rate matrix with n states
        AbstractRateMatrix(const AbstractRateMatrix& m);                                                                                //!< Copy constructor
        AbstractRateMatrix&                 operator=(const AbstractRateMatrix& r);                                                     //!< Assignment operator
        
        // protected methods available for derived classes
        std::vector<double>                 calculateStationaryFrequencies(void) const;                                                 //!< Calculate the stationary frequencies for the rate matrix
        bool                                checkTimeReversibity(double tolerance);
        virtual void                        computeStochasticMatrix(size_t n) const;
        virtual void                        computeDominatingRate(void) const;
        void                                exponentiateMatrixByScalingAndSquaring(double t,  TransitionProbabilityMatrix& p) const;
        
        // protected members available for derived classes
        MatrixReal*                         the_rate_matrix;                                                                            //!< Holds the rate matrix
        bool                                needs_update;
        
        // stochastic matrix
        mutable double                      dominating_rate;
        mutable std::vector<MatrixReal>     stochastic_matrix;                                                                          //!< Stochastic matrix raised to the power of n
    };
    
}

#endif

