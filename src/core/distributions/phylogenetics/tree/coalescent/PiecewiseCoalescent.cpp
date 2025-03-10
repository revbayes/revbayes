#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "PiecewiseCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "DemographicFunction.h"
#include "ConstantDemographicFunction.h"
#include "LinearDemographicFunction.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "RbMathLogic.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

/**
 * Default Constructor.
 *
 * @param N a vector a population sizes
 * @param i a vector of interval starts
 * @param meth the method for intervals. Options are 'EVENTS', 'SPECIFIED', 'UNIFORM'
 * @param c a vector of clade constraints
 *
 *
 *@note The parameter for interval starts, i, will not be used if the method for intervals is 'EVENTS' or 'UNIFORM'
 *
 *@note If the interval method is 'UNIFORM' then interval start times are equally distributed over the present time and the time of the root.
 *@note If the interval method is 'EVENTS' then we assume that the time between each coalescent event is an interval.
 *
 */
PiecewiseCoalescent::PiecewiseCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, const TypedDagNode<RbVector<long> > *n_events_pi, METHOD_TYPES meth, DEMOGRAPHY_FUNCTION_TYPES dem, const std::vector<Taxon> &tn, const std::vector<Clade> &c) :
    AbstractCoalescent( tn, c ),
    Nes( N ),
    interval_change_points_var( i ),
    number_events_per_interval( n_events_pi ),
    interval_method( meth ),
    demographic_function_var ( dem )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Nes );
    addParameter( interval_change_points_var );
    addParameter( number_events_per_interval );
    
    simulateTree();
    
    updateIntervals();
}



PiecewiseCoalescent::~PiecewiseCoalescent()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
PiecewiseCoalescent* PiecewiseCoalescent::clone( void ) const
{
    
    return new PiecewiseCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double PiecewiseCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double ln_prob = 0;
    
    // update the interval change points and coalescent sizes
    updateIntervals();
    
    // retrieved the coalescent ages
    std::vector<double> coalescent_ages;
    for (size_t i = 0; i < value->getNumberOfInteriorNodes()+1; ++i)
    {
        const TopologyNode& n = value->getInteriorNode( i );
        double a = n.getAge();
        coalescent_ages.push_back(a);
    }
    // sort the vector of ages in ascending order
    std::sort(coalescent_ages.begin(), coalescent_ages.end());
    
    // retrieve the ages of any serially sampled tips
    std::vector<double> serial_tip_ages;
    size_t num_taxa_at_present = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        const TopologyNode& n = value->getTipNode( i );
        double a = n.getAge();
        if ( a > 0.0 )
        {
            serial_tip_ages.push_back(a);
            --num_taxa_at_present;
        }
    }
    
    std::vector<double> combined_event_ages;
    std::vector<EVENT_TYPE> combined_event_types;

    // sort the vector of serial sampling times in ascending order
    std::sort(serial_tip_ages.begin(), serial_tip_ages.end());

    // initialize variables
    size_t at_coal_age                  = 0;
    size_t at_serial_age                = 0;
    size_t at_interval_change_point     = 0;
    double next_coal_age                = coalescent_ages[at_coal_age];
    double next_interval_change_point   = ( interval_change_points.size() > 0 ? interval_change_points[at_interval_change_point] : RbConstants::Double::inf );
    // set the next serially sampled tip age if present, otherwise Inf
    double next_serial_age              = (at_serial_age < serial_tip_ages.size() ? serial_tip_ages[at_serial_age] : RbConstants::Double::inf);

    // create master list of event times and types
    // events are either a sample (lineage size up), coalescence (lineage size down), or Ne changepoint (lineage size constant)
    do
    {
        next_coal_age = coalescent_ages[at_coal_age];
        if (next_coal_age <= next_serial_age && next_coal_age <= next_interval_change_point)
        {
            // coalescence
            combined_event_ages.push_back( next_coal_age );
            combined_event_types.push_back( COALESCENT );
            ++at_coal_age;
        }
        else if ( next_serial_age <= next_coal_age && next_serial_age <= next_interval_change_point )
        {
            // serial sample
            combined_event_ages.push_back( next_serial_age );
            combined_event_types.push_back( SERIAL_SAMPLE );
            ++at_serial_age;
            // set the next serial age
            if (at_serial_age < serial_tip_ages.size())
            {
                next_serial_age = serial_tip_ages[at_serial_age];
            }
            else
            {
                next_serial_age = RbConstants::Double::inf;
            }
        }
        else
        {
            // Ne change
            combined_event_ages.push_back( next_interval_change_point );
            combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
            ++at_interval_change_point;
            if ( at_interval_change_point < interval_change_points.size() )
            {
                next_interval_change_point = interval_change_points[at_interval_change_point];
            }
            else
            {
                next_interval_change_point = RbConstants::Double::inf;
            }
        }
    } while (at_coal_age < coalescent_ages.size());


    // size_t current_interval = 0;
    size_t current_num_lineages = num_taxa_at_present;
    double last_event_age = 0.0;

    size_t index_demographic_function = 0;

    for (size_t i = 0; i < combined_event_ages.size(); ++i)
    {
        
        double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
        double interval_area = getIntegral(last_event_age, combined_event_ages[i], index_demographic_function);
        
        // add log probability that nothing happens until the next event
        ln_prob -= n_pairs * interval_area;

        if (combined_event_types[i] == COALESCENT)
        {
            // coalescence
            double ne_at_coal_age = getDemographic(combined_event_ages[i], index_demographic_function);
            ln_prob -= log( ne_at_coal_age );
            --current_num_lineages;
        }
        else if (combined_event_types[i] == SERIAL_SAMPLE)
        {
            // sampled ancestor
            ++current_num_lineages;
        }
        else // combined_event_types[i] == DEMOGRAPHIC_MODEL_CHANGE
        {
            // Ne change, i.e., change of the demographic function
            ++index_demographic_function;
            if ( index_demographic_function > interval_change_points.size() )
            {
                throw RbException("Problem occurred in coalescent process with demographic functions: We tried to access a demographic function outside the vector.");
            }
        }

        last_event_age = combined_event_ages[i];
    }

    return ln_prob;
}


void PiecewiseCoalescent::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "getIntervalAges" )
    {
        
        updateIntervals();
        rv = interval_change_points;
        
    }
    else
    {
        throw RbException() << "The piecewise-constant coalescent process does not have a member method called '" << n << "'.";
    }
    
}


/**
 * Compute the population size of the demographic function at time event_age.
 */
double PiecewiseCoalescent::getDemographic(double event_age, double index) const
{
    
    // check which type of demographic function we have
    if ( demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::CONSTANT )
    {
        // get the constant demographic population size
        return pop_sizes[index];
    }
    else if (demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::LINEAR )
    {
        // get the linear demographic population size
        double N0 = pop_sizes[index];
        double N1 = ( index < (pop_sizes.size()-1) ? pop_sizes[index+1] : pop_sizes[index] );
        double t0 = ( index == 0 ? 0 : interval_change_points[index-1] );
        double t1 = ( index == interval_change_points.size() ? RbConstants::Double::inf : interval_change_points[index] );
        
        if ( t1 < t0 || t0 < 0 || N1 < 0 || event_age < t0 || event_age > t1)
        {
            throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
        }
        
        if ( N0 == N1 )
        {
            return N0;
        }
        else
        {
            double alpha = ( N1-N0 ) / (t1 - t0);
            return N0 + (event_age-t0) * alpha;
        }
    }
    else
    {
        throw RbException("Trying to compute demographic function value for not-implemented demographic function.");
    }
    
}


/**
 * Compute the integral over the population between the provided times and the given index of the demographic function.
 */
double PiecewiseCoalescent::getIntegral(double last_age, double event_age, double index) const
{
    
    // check which type of demographic function we have
    if ( demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::CONSTANT )
    {
        // get the constant demographic population size
        double delta = event_age - last_age;
        return delta / pop_sizes[index];
    }
    else if (demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::LINEAR )
    {
        // get the linear demographic population size
        double N0 = pop_sizes[index];
        double N1 = ( index < (pop_sizes.size()-1) ? pop_sizes[index+1] : pop_sizes[index] );
        double t0 = ( index == 0 ? 0 : interval_change_points[index-1] );
        double t1 = ( index == interval_change_points.size() ? RbConstants::Double::inf : interval_change_points[index] );
        
        if ( t1 < t0 || t0 < 0 || N1 < 0 || event_age < t0 || event_age > t1)
        {
            throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
        }
        
        if ( N0 == N1 )
        {
            double delta = event_age - last_age;
            return delta / N0;
        }
        else
        {
            double alpha = ( N1-N0 ) / (t1 - t0);
            return ( log( N0 + (event_age-t0) * alpha ) - log( N0 + (last_age-t0) * alpha ) ) / alpha;
        }
        
    }
    else
    {
        throw RbException("Trying to compute demographic function value for not-implemented demographic function.");
    }
}


/**
 * Compute the waiting time to the next coalescent event given the current age,
 * the random number (scaled by the current number of lineages), and
 * the index of the current interval/demographic function
 */
double PiecewiseCoalescent::getWaitingTime(double age, double rv, size_t index) const
{
    
    // check which type of demographic function we have
    if ( demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::CONSTANT )
    {
        // get the constant demographic population size
        return pop_sizes[index] * rv;

    }
    else if (demographic_function_var == DEMOGRAPHY_FUNCTION_TYPES::LINEAR )
    {
        // get the linear demographic population size
        double N0 = pop_sizes[index];
        double N1 = ( index < (pop_sizes.size()-1) ? pop_sizes[index+1] : pop_sizes[index] );
        double t0 = ( index == 0 ? 0 : interval_change_points[index-1] );
        double t1 = ( index == interval_change_points.size() ? RbConstants::Double::inf : interval_change_points[index] );
        
        if ( t1 < t0 || t0 < 0 || N1 < 0 || age < t0 || age > t1)
        {
            throw RbException("Impossible parameter values in Linear growth/decline demographic functions.");
        }
        
        if ( N0 == N1 )
        {
            return N0 * rv;
        }
        else
        {
            double alpha = ( N1-N0 ) / (t1 - t0);
            return (N0 + (age-t0) * alpha) * rv;
        }
        
    }
    else
    {
        throw RbException("Trying to compute demographic function value for not-implemented demographic function.");
    }
}




/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void PiecewiseCoalescent::keepSpecialization(const DagNode *affecter)
{
    
    // nothing to do here
    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PiecewiseCoalescent::restoreSpecialization(const DagNode *affecter)
{

    // Sebastian: This is currently redudant because we update the intervals each time when we compute the probability
    // just re-update the start times of the intervals
    updateIntervals();
    
}


/**
 * Simulate new coalescent times.
 *
 * \param[in]    n      The number of coalescent events to simulate.
 *
 * \return    A vector of the simulated coalescent times.
 */
std::vector<double> PiecewiseCoalescent::simulateCoalescentAges( size_t n ) const
{
    // first check if we want to set the interval ages here!
    if ( interval_method == SPECIFIED )
    {
        updateIntervals();
    }
    else
    {
        interval_change_points = RbVector<double>(Nes->getValue().size()-1, RbConstants::Double::inf);
        pop_sizes = Nes->getValue();
    }

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serial_tip_ages;
    size_t num_taxa_at_present = 0;
    for (size_t i = 0; i < num_taxa; ++i)
    {
        double a = taxa[i].getAge();
        if ( a > 0.0 )
        {
            serial_tip_ages.push_back(a);
        }
        else
        {
            ++num_taxa_at_present;
        }
    }

    
    // Put sampling times and pop-size changes into a single vector of event times
    std::vector<double>     combined_event_ages;
    std::vector<EVENT_TYPE> combined_event_types;
    
    // sort the vector of serial sampling times in ascending order
    std::sort(serial_tip_ages.begin(), serial_tip_ages.end());
    
    // initialize variables
    size_t at_coal_age                  = 0;
    size_t at_serial_age                = 0;
    size_t at_interval_change_point     = 0;
    double next_interval_change_point   = RbConstants::Double::inf;
    if ( interval_change_points.size() > 0 )
    {
        next_interval_change_point = interval_change_points[at_interval_change_point];
    }
    // set the next serially sampled tip age if present, otherwise Inf
    double next_serial_age              = (at_serial_age < serial_tip_ages.size() ? serial_tip_ages[at_serial_age] : RbConstants::Double::inf);

    
    // create master list of event times and types
    // pre-defined events are either a sample (lineage size up) or Ne changepoint (lineage size constant)
    size_t num_total_events = ( interval_method == SPECIFIED ? interval_change_points.size() + serial_tip_ages.size() : serial_tip_ages.size() );
    for (size_t num_events = 0; num_events < num_total_events; ++num_events)
    {
        if (next_interval_change_point <= next_serial_age)
        {
            // Ne change
            combined_event_ages.push_back( next_interval_change_point );
            combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE);
            ++at_interval_change_point;
            if (at_interval_change_point < interval_change_points.size())
            {
                next_interval_change_point = interval_change_points[at_interval_change_point];
            }
            else
            {
                next_interval_change_point = RbConstants::Double::inf;
            }
        }
        else
        {
            // serial sample
            combined_event_ages.push_back(next_serial_age);
            combined_event_types.push_back( SERIAL_SAMPLE );
            ++at_serial_age;
            if (at_serial_age < serial_tip_ages.size())
            {
                next_serial_age = serial_tip_ages[at_serial_age];
            }
            else
            {
                next_serial_age = RbConstants::Double::inf;
            }
        }
    }
 
    // cap vector with an event at t-infinity
    combined_event_ages.push_back( RbConstants::Double::inf );
    combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
    
    if ( interval_method == EVENTS )
    {
        // add dummies for the change point times
        for ( size_t i=0; i<interval_change_points.size(); ++i)
        {
            combined_event_ages.push_back( RbConstants::Double::inf );
            combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
        }
    }

    // allocate the vector for the times
    std::vector<double> coalescent_ages = std::vector<double>(n,0.0);
    
    size_t current_interval = 0;
    size_t ne_interval = 0;
    size_t current_num_lineages = num_taxa_at_present;
    
    size_t num_events_per_interval = size_t( ceil( double(n)/Nes->getValue().size()) );
    size_t current_num_events_in_interval = 0;

    
    // the current age of the process
    double sim_age = 0.0;

    // now simulate the ages
    for (size_t i = 0; i < n; ++i)
    {
        
        bool was_coalescent_event = false;
        do
        {
            double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
            double lambda = RbStatistics::Exponential::rv( n_pairs, *rng);
            double waitingTime = getWaitingTime(sim_age, lambda, ne_interval);
            
            if ( RbMath::isFinite(waitingTime) == false )
            {
                std::cerr << "Problem when computing waiting time." << std::endl;
            }
            
            sim_age += waitingTime;

            was_coalescent_event = (sim_age < combined_event_ages[current_interval] && current_num_lineages > 1) && waitingTime > 0;
            if ( was_coalescent_event == false )
            {
                // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time or Ne window, the number of active lineages changes
                // or the pop size changes, and it is necessary to discard any "excess" time,
                // which is drawn from an incorrect distribution,then we can draw a new time according to
                // the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current event in combinedEvents.
                sim_age = combined_event_ages[current_interval];
                if (combined_event_types[current_interval] == DEMOGRAPHIC_MODEL_CHANGE)
                {
                    // Ne change
                    ++ne_interval;
                }
                else
                {
                    // serial sample
                    ++current_num_lineages;
                }
                ++current_interval;
            }
        } while ( was_coalescent_event == false );

        coalescent_ages[i] = sim_age;
        
        ++current_num_events_in_interval;
        if ( interval_method == EVENTS && current_num_events_in_interval == num_events_per_interval && ne_interval < interval_change_points.size() )
        {
            interval_change_points[ne_interval] = sim_age;
            current_num_events_in_interval = 0;
            ++ne_interval;
        }
        
        --current_num_lineages;
    }
    
    return coalescent_ages;
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    old_p      Pointer to the old parameter.
 * \param[in]    new_p      Pointer to the new parameter.
 */
void PiecewiseCoalescent::swapParameterInternal(const DagNode *old_p, const DagNode *new_p)
{
    if (old_p == Nes)
    {
        Nes = static_cast<const TypedDagNode<RbVector<double> >* >( new_p );
    }
    else if (old_p == interval_change_points_var)
    {
        interval_change_points_var = static_cast<const TypedDagNode<RbVector<double> >* >( new_p );
    }
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PiecewiseCoalescent::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    // Sebastian: This is currently redudant because we update the intervals each time when we compute the probability
    // just update the start times of the intervals
//    updateIntervals();
    
}


/**
 * Recompute the current interval change point vector and corresponding population size vector.
 *
 *
 * @throw RbExpection when no interval start times are specified when the 'SPECIFIED' interval_method is used
 *
 */
void PiecewiseCoalescent::updateIntervals( void ) const
{
    
    // first check if we want to set the interval ages here!
    if ( interval_method == SPECIFIED )
    {
        // we obtained the interval ages, so we just need to make sure that they are sorted.
        if ( interval_change_points_var == NULL )
        {
            throw RbException("You have to provide the start times of the coalescent skyline intervals if you chose 'SPECIFIED' as the method of choice.");
        }
        
        // clean all the sets
        interval_change_points.clear();
        pop_sizes.clear();
                
        // important, we use this set of ages to automatically sort the ages.
        std::set<double> event_ages;
        
        const std::vector<double>& change_point_ages = interval_change_points_var->getValue();
        for (std::vector<double>::const_iterator it = change_point_ages.begin(); it != change_point_ages.end(); ++it)
        {
            event_ages.insert( *it );
        }
        
        // get our current population sizes
        const std::vector<double> &p = Nes->getValue();
                
        // we assume that the first population size corresponds to the interval starting at time 0
        pop_sizes.push_back( p[0] );
        
        size_t index = 0;
        size_t pos = 0;
        for (std::set<double>::iterator it = event_ages.begin(); it != event_ages.end(); ++it)
        {
            // get the age for the event
            double a = *it;
            
            // add the time to our vector
            interval_change_points.push_back( a );
            
            // add the population size at the event a
            pos = size_t( find(change_point_ages.begin(), change_point_ages.end(), a) - change_point_ages.begin() );
            if ( pos != change_point_ages.size() )
            {
                index = pos;
                pop_sizes.push_back( p[index+1] );
            }
            else
            {
                throw RbException("Couldn't sort vector of interval change points.");
            }
            
        }
        
        
    }
    else if ( interval_method == EVENTS )
    {
        // we directly use the population size given from the arguments
        pop_sizes = Nes->getValue();
        
        // next, we recompute the starting times of new intervals
        interval_change_points = RbVector<double>(Nes->getValue().size()-1, RbConstants::Double::inf);
        
        if ( this->value != NULL )
        {
            // retrieve the coalescent ages
            std::vector<double> coalescent_ages;
            for (size_t i = 0; i < value->getNumberOfInteriorNodes()+1; ++i)
            {
                const TopologyNode& n = value->getInteriorNode( i );
                double a = n.getAge();
                coalescent_ages.push_back(a);
            }
            // sort the vector of times in ascending order
            std::sort(coalescent_ages.begin(), coalescent_ages.end());
            
            // size_t num_events_per_interval = size_t( ceil( double(num_taxa-1.0)/Nes->getValue().size()) );
            size_t current_interval = 0;
            size_t current_num_events_in_interval = 0;
                        
            for (size_t i = 0; i < num_taxa-2; ++i)
            {
                ++current_num_events_in_interval;
                if ( current_num_events_in_interval == number_events_per_interval->getValue()[current_interval] )
                {
                    interval_change_points[current_interval] = coalescent_ages[i];
                    current_num_events_in_interval = 0;
                    ++current_interval;
                }
            }
        }

    }
    
}
