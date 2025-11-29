//
//  RateGenerator_Epoch.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__RateGenerator_Epoch__
#define __revbayes_proj__RateGenerator_Epoch__


#include "RateGenerator.h"
#include <complex>
#include <vector>
#include <map>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateGenerator_Epoch : public RateGenerator {
        
    public:
        RateGenerator_Epoch(size_t n, size_t ne);                                                                                                                   //!< Construct rate matrix with n states
//        RateGenerator_Epoch(const RateGenerator_Epoch& m);                                                                                                          //!< Copy constructor
        virtual                             ~RateGenerator_Epoch(void);                                                                                             //!< Destructor
        
        // overloaded operators
//        RateGenerator_Epoch&                operator=(const RateGenerator_Epoch& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateGenerator_Epoch*                clone(void) const;
        double                              getRate(size_t from, size_t to, double age, double rate) const;                                    //!< Calculate the rate from state i to state j over the given time interval scaled by a rate
        virtual RbVector<double>            getEpochTimesWithinInterval(double start_age, double end_age) const;
        const RbVector<RateGenerator>&      getRateGenerators(void) const;                                                                                         //!< Return the epoch generators
        const RbVector<double>&             getEpochTimes(void) const;                                                                                             //!< Return the epoch times
        const RbVector<double>&             getEpochRates(void) const;                                                                                             //!< Return the epoch rates
        void                                setEpochGenerators(const RbVector<RateGenerator>& rg);                                                              //!< Directly set the epoch generators
        void                                setEpochTimes(const RbVector<double> &t);                                                                           //!< Directly set the epoch times
        void                                setEpochRates(const RbVector<double>& r);                                                                           //!< Directly set the epoch rates
        virtual bool                        simulateStochasticMapping(double startAge, double endAge, double rate,std::vector<size_t>& transition_states, std::vector<double>& transition_times) const;
        virtual void                        printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const;            //!< print object for user (in user-formatted way)
        void                                update(void);
        
    private:
        size_t                              findEpochIndex( double t ) const;
        void                                assignEpochDominatingRates(void);
        void                                sampleBreakpointStates(std::vector<size_t>& breakpoint_states, std::vector<double> breakpoint_times, std::vector<TransitionProbabilityMatrix>& breakpoint_prob, double rate=1.0) const;
        void                                sampleNumberOfTransitionsPerInterval(std::vector<size_t>& num_events, std::vector<size_t> breakpoint_states, std::vector<double> breakpoint_times, std::vector<TransitionProbabilityMatrix> breakpoint_probs, std::vector<std::vector<MatrixReal> >& uniform_nth_power, double rate=1.0) const;

        
        // parameters
        RbVector<RateGenerator>             epochRateGenerators;
        RbVector<double>                    epochTimes;
        RbVector<double>                    epochRates;
        std::vector<double>                 epochDominatingRates;

        // helper variables
        size_t                              numEpochs;
        bool                                needs_update;
    };
    
    std::ostream&                       operator<<(std::ostream& o, const RateGenerator_Epoch& x);                                           //!< Overloaded output operator

}

#endif /* defined(__revbayes_proj__RateGenerator_Epoch__) */
