#include "DivergenceTimeCDF.h"

#include <cstdio>
#include <cmath>

#include "RbMathLogic.h"


using namespace RevBayesCore;


DivergenceTimeCDF::DivergenceTimeCDF(double o, double p, double u, double sf, const std::vector<double> &r, const std::vector<double> &b, const std::vector<double> &d) :
    origin_age( o ),
    present_age( p ),
    u( u ),
    rate_change_times( r ),
    birth( b ),
    death( d ),
    rho( sf )
{

}



double DivergenceTimeCDF::operator()(const double t)
{
    
    double p_0_T = 1.0 - pSurvival(present_age,origin_age,rho)  * exp( rateIntegral(present_age,origin_age) ) / rho;
    double p_0_t = 1.0 - pSurvival(origin_age-t,origin_age,rho) * exp( rateIntegral(origin_age-t,origin_age) ) / rho;
    double F_t   = p_0_t / p_0_T;
    
    return F_t - u;
}



double DivergenceTimeCDF::computeProbabilitySurvival(double start, double end) const
{
    // do the integration of int_{start}^{end} ( mu(s) exp(rate(t,s)) ds )
    // where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx ) - sum_{for all t < m_i < s in massExtinctionTimes }( log(massExtinctionSurvivalProbability[i]) )
    //
    // we compute the integral stepwise for each epoch
    
    double accummulated_rate_time = 0.0;
    double prev_time = start;
    double den = 1.0;
    
    size_t num_episodes = rate_change_times.size();
    double rate = 0.0;
    for ( size_t j=0; j<num_episodes; ++j )
    {
        // compute the rate
        rate = death[j] - birth[j];
        
        if ( (start < rate_change_times[j]) && (end >= rate_change_times[j]) )
        {
            // compute the integral for this time episode until the mass-extinction event
//            den += ( exp(-rate*prev_time) * death[j] / rate * exp( accummulated_rate_time ) * ( exp(rate* rate_change_times[j]) - exp(rate*prev_time)));
            den += ( death[j] / rate * exp( accummulated_rate_time ) * ( exp(rate * (rate_change_times[j]-prev_time)) - 1.0));

            if ( RbMath::isFinite(den) == false )
            {
                return 0.0;
            }
            
            accummulated_rate_time +=  (rate*(rate_change_times[j]-prev_time));
            // store the current time so that we remember from which episode we need to integrate next
            prev_time = rate_change_times[j];
            // integrate over the tiny time interval of the mass-extinction event itself and add it to the integral
            
        }
        
    }
    
    size_t index = 0;
    if ( rate_change_times.size() > 0 )
    {
        index = lower_index(end);
    }
    rate = death[index] - birth[index];
    
    // add the integral of the final epoch until the present time
    den = den + exp(-rate*prev_time) * exp( accummulated_rate_time ) * death[index] / rate * ( exp(rate*end) - exp(rate*prev_time));
    
    double res = 1.0 / den;
    
    return res;
}



/**
 *
 *
 */
size_t DivergenceTimeCDF::lower_index(double t) const
{
    return lower_index(t,0,rate_change_times.size());
}



/**
 *
 *
 */
size_t DivergenceTimeCDF::lower_index(double t, size_t min, size_t max) const
{
    size_t middle = min + (max - min) / 2;
    if ( min == middle )
    {
        return ( rate_change_times.size() > 0 && rate_change_times[max-1] < t ) ? max : min;
    }
    
    if ( rate_change_times[middle-1] > t )
    {
        return lower_index(t,min,middle);
    }
    else
    {
        return lower_index(t,middle,max);
    }
    
}

/**
 * Compute the probabililty of survival (no extinction) of the process including uniform taxon sampling at the present time.
 * The probability of survival is given by
 * [1 + int_{t_low}^{t_high} ( mu(s) exp(rate(t,s)) ds ) ]^{-1}
 * and can be simplified to
 * [1 + int_{t_low}^{t_high} ( mu'(s) exp(rate'(t,s)) ds ) - (r-1)/r*exp(rate'(t_low,t_high)) ]^{-1}
 * where mu' and rate' are the diversification rate function without incomplete taxon sampling.
 * Therefore we can just call pSurvival without incomplete taxon sampling that will be computed in the derived classes,
 * and add the sampling here so that sampling will be available for all models :)
 * For more information please read Hoehna, S. 2014. The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 * \param[in]    r          Sampling probability.
 *
 * \return The probability of survival of the process.
 */
double DivergenceTimeCDF::pSurvival(double start, double end, double r) const
{
    double rate = rateIntegral(start, end);
    double prob_surv = computeProbabilitySurvival(start, end);
    if ( prob_surv == 0.0 )
    {
        return 0.0;
    }
    else
    {
        double ps = 1.0 / prob_surv;
        
        return 1.0 / (ps - (r-1.0)/r * exp(rate) );
    }
}


double DivergenceTimeCDF::rateIntegral(double start, double end) const
{
    
    double accummulated_rate_time = 0.0;
    double prev_time = start;
    
    size_t num_episodes = rate_change_times.size();
    double rate = 0.0;
    for ( size_t j=0; j<num_episodes; ++j )
    {
        // compute the rate
        rate = death[j] - birth[j];
        
        if ( (start < rate_change_times[j]) && (end >= rate_change_times[j]) )
        {
            // compute the integral for this time episode until the mass-extinction event
            
            accummulated_rate_time +=  (rate*(rate_change_times[j]-prev_time));
            // store the current time so that we remember from which episode we need to integrate next
            prev_time = rate_change_times[j];
            // integrate over the tiny time interval of the mass-extinction event itself and add it to the integral
            
        }
        
    }
    
    size_t index = 0;
    if ( rate_change_times.size() > 0 )
    {
        index = lower_index(end);
    }
    rate = death[index] - birth[index];
    accummulated_rate_time += (rate*(end-prev_time));
    
    return accummulated_rate_time;
}

