#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

#include "Clade.h"
#include "DivergenceTimeCDF.h"
#include "EpisodicBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathLogic.h"
#include "BirthDeathProcess.h"
#include "RbException.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"

#include "boost/format.hpp" // IWYU pragma: keep
#include "boost/math/tools/toms748_solve.hpp"
#include "boost/optional/optional.hpp"
#include <boost/math/tools/roots.hpp> // IWYU pragma: keep


namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

EpisodicBirthDeathProcess::EpisodicBirthDeathProcess(const TypedDagNode<double> *ra,
                                                     const TypedDagNode<RbVector<double> > *sr,
                                                     const TypedDagNode<RbVector<double> > *st,
                                                     const TypedDagNode<RbVector<double> > *er,
                                                     const TypedDagNode<RbVector<double> > *et,
                                                     const TypedDagNode<double> *r,
                                                     const TypedDagNode<double> *mp,
                                                     const std::string& ss,
                                                     const std::vector<Clade> &ic,
                                                     const std::string &cdt,
                                                     const std::vector<Taxon> &tn,
                                                     Tree* t) : BirthDeathProcess( ra, r, mp, ss, ic, cdt, tn, t ),
    lambda_rates( sr ),
    lambda_times( st ),
    mu_rates( er ),
    mu_times( et )
{
    addParameter( lambda_rates );
    addParameter( lambda_times );
    addParameter( mu_rates );
    addParameter( mu_times );
    
    prepareProbComputation();

    if ( starting_tree == NULL )
    {
        simulateTree(true);
    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
EpisodicBirthDeathProcess* EpisodicBirthDeathProcess::clone( void ) const
{
    
    return new EpisodicBirthDeathProcess( *this );
}


double EpisodicBirthDeathProcess::computeProbabilitySurvival(double start, double end) const
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
//            den += ( exp(-rate*prev_time) * death[j] / rate * exp( accummulated_rate_time ) * ( exp(rate*rate_change_times[j]) - exp(rate*prev_time)));
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



double EpisodicBirthDeathProcess::lnSpeciationRate(double t) const
{
    size_t index = lower_index(t);
    double ln_lambda = log( birth[ index ] );
    
    return ln_lambda;
}


/**
 *
 *
 */
size_t EpisodicBirthDeathProcess::lower_index(double t) const
{
    return lower_index(t,0,rate_change_times.size());
}



/**
 *
 *
 */
size_t EpisodicBirthDeathProcess::lower_index(double t, size_t min, size_t max) const
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
 *
 *
 */
void EpisodicBirthDeathProcess::prepareProbComputation( void ) const
{
    
    // clean all the sets
    rate_change_times.clear();
    birth.clear();
    death.clear();
    
    double present_time = process_age->getValue();
    
    std::set<double> event_times;
    
    const std::vector<double>& birth_times = lambda_times->getValue();
    for (std::vector<double>::const_iterator it = birth_times.begin(); it != birth_times.end(); ++it)
    {
        event_times.insert( *it );
    }
    
    const std::vector<double>& death_times = mu_times->getValue();
    for (std::vector<double>::const_iterator it = death_times.begin(); it != death_times.end(); ++it)
    {
        event_times.insert( *it );
    }
    
    const std::vector<double> &b = lambda_rates->getValue();
    const std::vector<double> &d = mu_rates->getValue();
    
    size_t index_birth = b.size()-1;
    size_t index_death = d.size()-1;
    
    birth.push_back( b[index_birth] );
    death.push_back( d[index_death] );
    
    size_t pos = 0;
    for (std::set<double>::reverse_iterator it = event_times.rbegin(); it != event_times.rend(); ++it)
    {
        double t = *it;
        
        // add the time to our vector
        rate_change_times.push_back( present_time - t );
        
        // add the speciation rate at the rate-change event t
        pos = size_t( find(birth_times.begin(), birth_times.end(), t) - birth_times.begin() );
        if ( pos != birth_times.size() )
        {
            index_birth = pos;
            birth.push_back( b[index_birth] );
        }
        else
        {
            birth.push_back( b[index_birth] );
        }
        
        
        // add the extinction rate at the rate-change event t
        pos = size_t( find(death_times.begin(), death_times.end(), t) - death_times.begin() );
        if ( pos != death_times.size() )
        {
            index_death = pos;
            death.push_back( d[index_death] );
        }
        else
        {
            death.push_back( d[index_death] );
        }
        
    }
    
}


void EpisodicBirthDeathProcess::prepareSurvivalProbability(double end, double r)
{
    // do the integration of int_{start}^{end} ( mu(s) exp(rate(t,s)) ds )
    // where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx ) - sum_{for all t < m_i < s in massExtinctionTimes }( log(massExtinctionSurvivalProbability[i]) )
    //
    // we compute the integral stepwise for each epoch
    
    
    double accummulated_rate_time_until_present = 0.0;
    double den = 0.0;
    
    
    // get the number of epochs
    size_t num_epochs = rate_change_times.size();
    
    int j=int(num_epochs)-1;
    while ( j >= 0 && end < rate_change_times[j] )
    {
        --j;
    }
    
    // compute the rate
    double rate = death[j+1] - birth[j+1];
    
    size_t num_taxa = value->getNumberOfTips();

    log_p_survival  = std::vector<double>(num_taxa-2,0.0);

    double prev_end = end;
    for ( size_t i=num_taxa-2; i>0; --i )
    {
        double start = divergence_times[i];
        
        double factor = 1.0;
        double summand = 0.0;
    
        while ( j >= 0 )
        {
            
            if ( start < rate_change_times[j] )
            {
                // compute the integral for this time episode until the rate-shift event
                double delta = prev_end-rate_change_times[j];
                factor = exp(rate*delta);
                den *= factor;
                
                accummulated_rate_time_until_present += rate * delta;

                summand = ( exp(-rate*rate_change_times[j]) * death[j+1] / rate * ( exp(rate*prev_end) - exp(rate*rate_change_times[j]) ) );
                den += summand;
                
                // store the current time so that we remember from which epoch we need to integrate next
                prev_end = rate_change_times[j];

                --j;
                
                // re-compute the rate
                rate = death[j+1] - birth[j+1];
                
            }
            else
            {
                break;
            }
            
        }
        
        // rescale the previous sum of exponentiated rates
        factor = exp(rate*(prev_end-start));
        den *= factor;

        summand = exp(-rate*start) * death[j+1] / rate * ( exp(rate*prev_end) - exp(rate*start) );
        den += summand;

        double ps = 1.0 + den;
        
        accummulated_rate_time_until_present += rate * (prev_end-start);
        log_p_survival[i-1] = -log(ps - (r-1.0)/r * exp(accummulated_rate_time_until_present) );

        prev_end = start;
        
    }
    
}


void EpisodicBirthDeathProcess::prepareRateIntegral(double end)
{
    
    double accummulated_rate_time = 0.0;
    
    size_t num_episodes = rate_change_times.size();
    double rate = 0.0;
    
    int j=int(num_episodes)-1;
    while ( j >= 0 && end < rate_change_times[j] )
    {
        --j;
    }
    
    // compute the rate
    rate = death[j+1] - birth[j+1];
    
    size_t num_taxa = value->getNumberOfTips();

    rate_integral  = std::vector<double>(num_taxa-2,0.0);

    double prev_end = end;
    for ( size_t i=num_taxa-2; i>0; --i )
    {
        double start = divergence_times[i];
        
        while ( j >= 0 )
        {
        
            if ( start < rate_change_times[j] )
            {
                // compute the integral for this time episode until the rate-shift event
                accummulated_rate_time += (rate*(prev_end-rate_change_times[j]));
                // store the current time so that we remember from which episode we need to integrate next
                prev_end = rate_change_times[j];
                
                --j;
                
                // re-compute the rate
                rate = death[j+1] - birth[j+1];
            }
            else
            {
                break;
            }
            
        }
        
        accummulated_rate_time += (rate*(prev_end-start));
        rate_integral[i-1] = accummulated_rate_time;
        
        prev_end = start;
        
    }
    
}


double EpisodicBirthDeathProcess::rateIntegral(double start, double end) const
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



double EpisodicBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant EBDP
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    size_t index = 0;
    if ( rate_change_times.size() > 0 )
    {
        index = lower_index(present);
    }

    // get the parameters
    double age = origin - present;
    double b = birth[index];
    double d = death[index];
    double r = rho->getValue();
    
    
    // get a random draw
    double u = rng->uniform01();
    
    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    else
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    

    DivergenceTimeCDF cdf = DivergenceTimeCDF(age, present, u, r, rate_change_times, birth, death);
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=100;
    boost::math::tools::eps_tolerance<double> tol(10);
//    DivergenceTimeCDF_tolerance<double> tol( (age - present) / 1000.0);
    
    double low = cdf.operator()(present);
    double high = cdf.operator()(age);
    
    if ( RbMath::isFinite( low ) == false || RbMath::isFinite(high) == false )
    {
        
        throw RbException("Cannot simulate a tree under the EBD model for the given parameters.");
        
        double p_0_T = 1.0 - pSurvival(present,age,r) * exp( rateIntegral(present,age) ) / r;
        double p_0_t = 1.0 - pSurvival(age-t,age,r) * exp( rateIntegral(age-t,age) ) / r;
        double F_t   = p_0_t / p_0_T;
        
        double a1 = pSurvival(present,age,r);
        double a2 = exp( rateIntegral(present,age) );
        double a3 = pSurvival(age-t,age,r);
        double a4 = exp( rateIntegral(age-t,age) );
        double a5 = exp( lnProbSurvival(present,age,r) + rateIntegral(present, age) );
        double a6 = exp( lnProbSurvival(age-t,age,r) + rateIntegral(age-t, age) );
        
        t = present;
        
        double p_0_T_low = 1.0 - pSurvival(present,age,r) * exp( rateIntegral(present,age) ) / r;
        double p_0_t_low = 1.0 - pSurvival(age-t,age,r) * exp( rateIntegral(age-t,age) ) / r;
        double F_t_low   = p_0_t / p_0_T;
        
        t = age;
        
        double p_0_T_high = 1.0 - pSurvival(present,age,r) * exp( rateIntegral(present,age) ) / r;
        double p_0_t_high = 1.0 - pSurvival(age-t,age,r) * exp( rateIntegral(age-t,age) ) / r;
        double F_t_high   = p_0_t / p_0_T;
    }
    
    Result r1 = boost::math::tools::toms748_solve(cdf, present, age, tol, max_iter);
//    Result r1 = boost::math::tools::bisect(cdf, present, age, tol, max_iter);

    return present + (r1.first+r1.second)/2.0;
}


/** Swap a parameter of the distribution */
void EpisodicBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == lambda_rates )
    {
        lambda_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == lambda_times )
    {
        lambda_times = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == mu_rates )
    {
        mu_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == mu_times )
    {
        mu_times = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else
    {
        // delegate the super-class
        BirthDeathProcess::swapParameterInternal(oldP, newP);
    }
    
}
