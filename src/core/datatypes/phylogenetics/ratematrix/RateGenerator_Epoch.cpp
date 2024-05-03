//
//  RateGenerator_Epoch.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include <iomanip>
#include <cstddef>
#include <iostream>
#include <vector>

#include "AbstractRateMatrix.h"
#include "DistributionPoisson.h"
#include "RateGenerator_Epoch.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "RateMatrix.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateGenerator_Epoch::RateGenerator_Epoch(size_t n, size_t ne) : RateGenerator( n ), numEpochs( ne ), needs_update(true)
{
    epochTimes = RbVector<double>(1, 0.0);
    epochRates = RbVector<double>(1, 1.0);
    epochRateGenerators = RbVector<RateGenerator>(1,  RateMatrix_JC(this->num_states));
    update();
    
}

/** Destructor */
RateGenerator_Epoch::~RateGenerator_Epoch(void)
{

}

/** Calculate the transition probabilities */
void RateGenerator_Epoch::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    // what amount of tole
    double precisionError = 1E-6;
    double diffAge = startAge - endAge;
    
    if (diffAge + precisionError < 0)
    {
        throw RbException("RateGenerator_Epoch cannot compute transition probabilities for negative branch lengths");
    }
    
    // P = I
    TransitionProbabilityMatrix tp(num_states);
    for (size_t i = 0; i < num_states; i++)
        tp[i][i] = 1.0;
    
    if (diffAge > 0)
    {
        // get current time on branch
        double currAge = startAge;
        
        // find epoch
        size_t epochIdx = findEpochIndex(currAge);
        
        // multiply transition probs across epochs
        while (currAge > endAge)
        {
            double nextAge = 0.0;
            if (epochIdx < numEpochs)
                nextAge = epochTimes[epochIdx];
            
            // get next time, which is the next epoch or branch end
            if (nextAge < endAge)
                nextAge = endAge;
    //        double dt = currAge - nextAge;
            
            // first, get the rate matrix for this branch
            const RateGenerator& rg = epochRateGenerators[epochIdx];
            double r = epochRates[epochIdx];
            rg.calculateTransitionProbabilities( currAge, nextAge, r * rate, P );

            // epochs construct DTMC
            tp *= P;
            
//            std::cout << P << "\n\n";
//            std::cout << tp << "\n";
//            std::cout << "-------\n";
            
            // advance increment
            currAge = nextAge;
            epochIdx++;
        }
    }
    else
    {
        ; // do nothing, i.e. transition probabilty matrix equals the identity matrix
    }
    
    P = TransitionProbabilityMatrix(tp);
//    std::cout << P << "\n";
}


RateGenerator_Epoch* RateGenerator_Epoch::clone( void ) const
{
    
    return new RateGenerator_Epoch( *this );
}


size_t RateGenerator_Epoch::findEpochIndex( double t ) const
{
    int i;
   
    for (i = 0; i < epochTimes.size() && t < epochTimes[i]; i++)
    {
        ; // do nothing
    }
    return i;
}


double RateGenerator_Epoch::getRate(size_t from, size_t to, double age, double rate) const
{
    size_t epochIdx = findEpochIndex(age);
   
    const RateGenerator *rg = &epochRateGenerators[epochIdx];

    double r = rg->getRate(from, to, age, rate * epochRates[epochIdx]);
    
    return r;
}

RbVector<double> RateGenerator_Epoch::getEpochTimesWithinInterval(double start_age, double end_age) const
{
    RbVector<double> v;
    size_t epoch_idx = findEpochIndex(start_age);
    while( epochTimes[epoch_idx] > end_age ) {
        v.push_back( epochTimes[epoch_idx] );
        epoch_idx++;
    }
    v.push_back(end_age);
    
    return v;
}

const RbVector<RateGenerator>& RateGenerator_Epoch::getRateGenerators(void) const
{
    return epochRateGenerators;
}

const RbVector<double>& RateGenerator_Epoch::getEpochTimes(void) const
{
    return epochTimes;
}

const RbVector<double>& RateGenerator_Epoch::getEpochRates(void) const
{
    return epochRates;
}


void RateGenerator_Epoch::assignEpochDominatingRates(void)
{
    epochDominatingRates = std::vector<double>(epochTimes.size());
    for (size_t i = 0; i < epochTimes.size(); i++)
    {
        double max_rate = 0.0;
        for (size_t j = 0; j < epochRateGenerators[0].size(); j++)
        {
            double r = -1.0 * getRate(j, j, epochTimes[i], 1.0);
            if (r > max_rate) {
                max_rate = r;
            }
        }
        
        epochDominatingRates[i] = max_rate;
    }
    
}

void RateGenerator_Epoch::setEpochGenerators(const RbVector<RateGenerator>& rg)
{
    epochRateGenerators = rg;
    needs_update = true;
    assignEpochDominatingRates();
}

void RateGenerator_Epoch::setEpochTimes(const RbVector<double> &t)
{
    epochTimes = t;
    needs_update = true;
}

void RateGenerator_Epoch::setEpochRates(const RbVector<double>& r)
{
    epochRates = r;
    needs_update = true;
    assignEpochDominatingRates();
}

bool RateGenerator_Epoch::simulateStochasticMapping(double startAge, double endAge, double rate, std::vector<size_t>& transition_states, std::vector<double>& transition_times) const
{
    if ( dynamic_cast<const RateMatrix*>( &epochRateGenerators[0] ) == NULL )
    {
        throw RbException("RateGenerator_Epoch::simulateStochasticMapping only supported when epoch rates are defined by RateMatrix objects");
    }
    
    bool success = false;

    // Trees may have branches of zero length, e.g. sampled ancestors under FBDP
    // These branches always have start states that match end states.
    if (startAge == endAge)
    {
        transition_times.push_back(0.0);
        return true;
    }

    // compute the breakpoint times
    std::vector<double> breakpoint_times;
    std::vector<size_t> breakpoint_index;
    breakpoint_times.push_back(startAge);
    breakpoint_index.push_back( findEpochIndex(startAge) );
    for (size_t i = 0; i < epochTimes.size(); i++) {
        if (epochTimes[i] < startAge && epochTimes[i] >= endAge) {
            breakpoint_times.push_back(epochTimes[i]);
            breakpoint_index.push_back(i+1);
        }
        else if (epochTimes[i] < endAge)
            break;
    }
    if (breakpoint_times.back() != endAge) {
        breakpoint_times.push_back(endAge);
    }
    
    // sample epoch breakpoint states
    std::vector<size_t> breakpoint_states = transition_states;
    std::vector<TransitionProbabilityMatrix> breakpoint_probs;
    sampleBreakpointStates(breakpoint_states, breakpoint_times, breakpoint_probs, rate);

    // make stochastic map for each interval
    transition_states.clear();
    transition_times.clear();
    for (size_t i = 1; i < breakpoint_times.size(); i++)
    {
        // interval variables
        size_t interval_index = breakpoint_index[i-1];
        double interval_start_age = breakpoint_times[i-1];
        double interval_end_age = breakpoint_times[i];
        double interval_rate = rate * epochRates[interval_index];
        
        std::vector<size_t> interval_transition_states;
        interval_transition_states.push_back(breakpoint_states[i-1]);
        interval_transition_states.push_back(breakpoint_states[i]);
        std::vector<double> interval_transition_times;
        
        // generate stochastic mapping for interval
        const RateMatrix* rate_matrix = dynamic_cast<const RateMatrix*>( &epochRateGenerators[interval_index] );
        success |= const_cast<RateMatrix*>(rate_matrix)->simulateStochasticMapping(interval_start_age, interval_end_age, interval_rate, interval_transition_states, interval_transition_times);
        
        // add the new events
        transition_states.insert( transition_states.end(), interval_transition_states.begin(), interval_transition_states.end() );
        transition_times.insert( transition_times.end(), interval_transition_times.begin(), interval_transition_times.end() );
    }
    
    // filter out the virtual events
    std::vector<size_t> save_states(1, transition_states[0]);
    std::vector<double> save_times(1, transition_times[0]);
    size_t prev_state = transition_states[0];
    for (size_t i = 1; i < transition_states.size(); i++) {
        if (prev_state != transition_states[i]) {
            // new real event detected
            save_states.push_back(transition_states[i]);
            save_times.push_back(transition_times[i]);
            prev_state = transition_states[i];
        }
        else
        {
            // add duration for virtual event to last real event
            size_t prev_idx = save_states.size() - 1;
            save_times[prev_idx] += transition_times[i];
        }
    }
    
    transition_states = save_states;
    transition_times = save_times;

    return success;
}

void RateGenerator_Epoch::sampleBreakpointStates(std::vector<size_t>& breakpoint_states, std::vector<double> breakpoint_times, std::vector<TransitionProbabilityMatrix>& breakpoint_probs, double rate) const
{
    size_t num_breaks = breakpoint_times.size() - 2;
    size_t num_intervals = num_breaks + 1;
    
    
    
    // transition probabilities for the interval (t_{i-1}, t_{i})
    for (size_t i = 1; i < breakpoint_times.size(); i++)
    {
        TransitionProbabilityMatrix P(this->num_states);
        calculateTransitionProbabilities(breakpoint_times[i-1], breakpoint_times[i], rate, P);
        breakpoint_probs.push_back(P);
    }
    
    // transition probabilities for the interval (t_{i}, t_{n})
    std::vector<TransitionProbabilityMatrix> prob_n;
    for (size_t i = 0; i < num_intervals; i++)
    {
        TransitionProbabilityMatrix P(this->num_states);
        calculateTransitionProbabilities(breakpoint_times[i], breakpoint_times[num_intervals], rate, P);
        prob_n.push_back(P);
    }
    
    // sample states per epoch breakpoint
    size_t end_state = breakpoint_states[1];
    for (size_t i = 1; i < num_intervals; i++)
    {
        // find the state that precedes interval i's state
        size_t prev_state = breakpoint_states[ breakpoint_states.size()-2 ];
        
        // compute probs of sampling intermediate state s
        std::vector<double> prob_curr_state(this->num_states, 0.0);
        double sum_curr_state = 0.0;
        for (size_t s = 0; s < this->num_states; s++)
        {
            double p_is = breakpoint_probs[i-1][prev_state][s];
            double p_sn = prob_n[i][s][end_state];
            double p_in = prob_n[i-1][prev_state][end_state];
            // std::cout << prev_state << " -> " << s << " -> " << end_state << " " << p_is << " " << p_sn << " " << p_in << "\n";
            prob_curr_state[s] = p_is * p_sn / p_in;
            sum_curr_state += prob_curr_state[s];
        }
        
        // sample breakpoint state
        double u = GLOBAL_RNG->uniform01();
        size_t curr_state = 0;
        while (u > 0.0) {
            u -= prob_curr_state[curr_state];
            if (u <= 0.0) {
                break;
            }
            curr_state += 1;
        }
        
        // insert the sampled state
        breakpoint_states.insert( breakpoint_states.end()-1, curr_state );
    }
}


void RateGenerator_Epoch::sampleNumberOfTransitionsPerInterval(std::vector<size_t> &num_events,
                                                               std::vector<size_t> breakpoint_states,
                                                               std::vector<double> breakpoint_times,
                                                               std::vector<TransitionProbabilityMatrix> breakpoint_probs,
                                                               std::vector<std::vector<MatrixReal> >& uniform_nth_power,
                                                               double rate) const
{
    
    size_t num_segments = breakpoint_times.size() - 1;
    
    // set up uniformized matrices
    std::vector<MatrixReal> uniform_first_power;
//    std::vector<std::vector<MatrixReal> > uniform_nth_power( num_segments );
    
    // build the base uniformization objects
    MatrixReal identity_matrix(this->num_states,this->num_states);
    for (size_t i = 0; i < this->num_states; i++)
        identity_matrix[i][i] = 1.0;
    
    for (size_t i = 0; i < num_segments; i++)
    {
        size_t epoch_idx = findEpochIndex(breakpoint_times[i+1]);
        const AbstractRateMatrix* rm = dynamic_cast<const AbstractRateMatrix*>( &epochRateGenerators[epoch_idx] );
        double r = epochDominatingRates[epoch_idx];
        uniform_first_power.push_back( rm->getRateMatrix() );
        uniform_first_power[i] *= (1.0 / r );
        for (size_t j = 0; j < this->num_states; j++)
        {
            uniform_first_power[i][j][j] += 1.0;
        }
    }
    
    // sample from f(n|a,b,lambda)
    for (size_t i = 0; i < num_segments; i++)
    {
        // segment length and start/end states
        size_t epoch_idx = findEpochIndex(breakpoint_times[i]);
        double dt = (breakpoint_times[i] - breakpoint_times[i+1]);
        double r = rate * epochRates[epoch_idx] * epochDominatingRates[epoch_idx];
        size_t a = breakpoint_states[i];
        size_t b = breakpoint_states[i+1];
        
        
//        std::cout << a << " -> " << b << "\n";
//        std::cout << breakpoint_times[i] << " " << breakpoint_times[i+1] << "\n";
//        std::cout << dt << "\n";
//        std::cout << r << "\n";
//        std::cout << epoch_idx << "\n";
        
        // sample cmf
        double u = GLOBAL_RNG->uniform01();
        double p_ab = breakpoint_probs[i][a][b];
        double g = p_ab * u;
        
        size_t n = 0;
        while (g > 0.0) {
            
            // prob of sampling n events
            double p_n = RbStatistics::Poisson::pdf(r*dt, (int)n);
            
            if (n == 0) {
                uniform_nth_power[i].push_back( identity_matrix );
            } else if (n == 1) {
                uniform_nth_power[i].push_back( uniform_first_power[i] );
            } else {
                uniform_nth_power[i].push_back( uniform_nth_power[i][n-1] * uniform_first_power[i] );
            }
//            std::cout << uniform_nth_power[i][n] << "\n";
            
            // prob of a->b after n events in segment i
            double p_b_given_a = uniform_nth_power[i][n][a][b];
            
            // unnormalized total prob
            double p_total = p_n * p_b_given_a;
            
            // sample n events
            g -= p_total;
            if (g <= 0.0) {
                break;
            }
            n += 1;
            
            if (n > 100) {
                
                throw RbException("RateGenerator_Epoch::sampleNumberOfTransitionsPerInterval failed when trying to sample >100 events!");
            }
        }
        num_events.push_back(n);
    }
}

void RateGenerator_Epoch::update( void ) {
    
    if ( needs_update )
    {
        // clean flags
        needs_update = false;
    }
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RateGenerator_Epoch& x) {
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << "[ ";
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < x.size(); i++) {
        if (i == 0)
            o << "[ ";
        else
            o << "  ";
        
        for (size_t j = 0; j < x.size(); ++j)
        {
            if (j != 0)
                o << ", ";
        }
        o <<  " ]";
        
        if (i == x.size()-1)
            o << " ]";
        else
            o << " ,\n";
        
    }
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
    return o;
}

void RateGenerator_Epoch::printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const {
    
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    std::stringstream prev_time( "Inf" );
    for (size_t i = 0; i < epochTimes.size(); i++) {
        
        std::stringstream curr_time;
        curr_time << epochTimes[i];
        o << "\n";
        o << "Epoch with index " << i << "\n";
        o << "Epoch time from " << prev_time.str() << " to " << curr_time.str() << "\n";
        o << std::setprecision(4) << std::fixed;
        o << "Epoch rate multiplier " << epochRates[i] << "\n";
        o.setf(previousFlags);
        o.precision(previousPrecision);
        o << "Epoch rate generator:\n";
        // append each clado matrix's print statement
        epochRateGenerators[i].printForUser(o, sep, l, left);
        o << "\n";
        prev_time.str( curr_time.str() );
    }
}
