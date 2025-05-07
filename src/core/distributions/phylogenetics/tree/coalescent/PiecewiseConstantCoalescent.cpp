#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "PiecewiseConstantCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

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
PiecewiseConstantCoalescent::PiecewiseConstantCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, METHOD_TYPES meth, const std::vector<Taxon> &tn, const std::vector<Clade> &c) :
    AbstractCoalescent( tn, c ),
    Nes( N ),
    interval_starts_var( i ),
    interval_method( meth )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Nes );
    addParameter( interval_starts_var );
    
    simulateTree();
    
    updateIntervals();
}



PiecewiseConstantCoalescent::~PiecewiseConstantCoalescent()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
PiecewiseConstantCoalescent* PiecewiseConstantCoalescent::clone( void ) const
{
    
    return new PiecewiseConstantCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double PiecewiseConstantCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    
    const_cast< PiecewiseConstantCoalescent* >( this )->updateIntervals();
    
    // retrieved the speciation times
    std::vector<double> ages;
    for (size_t i = 0; i < value->getNumberOfInteriorNodes()+1; ++i)
    {
        const TopologyNode& n = value->getInteriorNode( i );
        double a = n.getAge();
        ages.push_back(a);
    }
    // sort the vector of times in ascending order
    std::sort(ages.begin(), ages.end());
    
    size_t currentInterval = 0;
    
    for (size_t i = 0; i < ages.size(); ++i)
    {
        size_t j = num_taxa - i;
        double theta = 0.0;
        double nPairs = j * (j-1) / 2.0;
        
        double prevCoalescentTime = 0.0;
        if ( i > 0 )
        {
            prevCoalescentTime = ages[i-1];
        }
        
        double deltaAge = ages[i];
        bool valid = false;
        do
        {
            theta = pop_sizes[currentInterval];
            double age = ages[i];
            double max = age;
            if ( currentInterval < interval_starts.size() )
            {
                max = (age > interval_starts[currentInterval]) ? interval_starts[currentInterval] : age;
            }
            
            deltaAge = max - prevCoalescentTime;
            valid = currentInterval >= interval_starts.size() || age <= interval_starts[currentInterval];
            if ( !valid )
            {
                ++currentInterval;
                prevCoalescentTime = max;
                
                lnProbTimes -= nPairs * deltaAge / theta ;
            }
            
        } while ( !valid );
        
        lnProbTimes += log( 1.0 / theta ) - nPairs * deltaAge / theta;
    }
    
    return lnProbTimes;
}


void PiecewiseConstantCoalescent::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "getIntervalAges" )
    {
        
        rv = interval_starts;
        
    }
//    else if ( n == "averageExtinctionRate" )
//    {
//
//
//    }
    else
    {
        throw RbException("The piecewise-constant coalescent process does not have a member method called '" + n + "'.");
    }
    
}


/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void PiecewiseConstantCoalescent::keepSpecialization(const DagNode *affecter)
{
    
    // nothing to do here
    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PiecewiseConstantCoalescent::restoreSpecialization(const DagNode *affecter)
{
    
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
std::vector<double> PiecewiseConstantCoalescent::simulateCoalescentAges( size_t n ) const
{
    
    // first check if we want to set the interval ages here!
    if ( interval_method == SPECIFIED )
    {
        if ( interval_starts_var != NULL )
        {
            interval_starts = interval_starts_var->getValue();
        }
        else
        {
            throw RbException("You have to provide the start times of the coalescent skyline intervals if you chose 'SPECIFIED' as the method of choece.");
        }
        const_cast< PiecewiseConstantCoalescent* >( this )->updateIntervals();
    }
    else
    {
        interval_starts = RbVector<double>(Nes->getValue().size()-1, RbConstants::Double::inf);
    }
    size_t num_events_per_interval = size_t( ceil( double(n)/Nes->getValue().size()) );
    size_t current_num_events_in_interval = 0;
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // now simulate the ages
    
    // allocate the vector for the times
    std::vector<double> coalescent_times = std::vector<double>(n,0.0);
    
    
    const RbVector<double> &pop_sizes  = Nes->getValue();
    size_t current_interval = 0;
    // draw a time for each speciation event condition on the time of the process
    for (size_t i = 0; i < n; ++i)
    {
        double prevCoalescentTime = 0.0;
        if ( i > 0 )
        {
            prevCoalescentTime = coalescent_times[i-1];
        }
        
        size_t j = num_taxa - i;
        double nPairs = j * (j-1) / 2.0;
        
        double sim_age = 0.0;
        bool valid = false;
        do
        {
//            double theta = 1.0 / (2.0*popSizes[current_interval]);
//            double lambda = nPairs * theta / 2.0;
            double theta = 1.0 / (pop_sizes[current_interval]);
            double lambda = nPairs * theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            sim_age = prevCoalescentTime + u;
            valid = current_interval >= interval_starts.size() || sim_age < interval_starts[current_interval];
            if ( valid == false )
            {
                prevCoalescentTime = interval_starts[current_interval];
                ++current_interval;
            }
            
        } while ( valid == false );
        
        coalescent_times[i] = sim_age;
        
        ++current_num_events_in_interval;
        if ( interval_method == EVENTS && current_num_events_in_interval == num_events_per_interval && current_interval < interval_starts.size() )
        {
            interval_starts[current_interval] = sim_age;
            current_num_events_in_interval = 0;
        }
        
    }
    
    return coalescent_times;
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void PiecewiseConstantCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == Nes)
    {
        Nes = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    else if (oldP == interval_starts_var)
    {
        interval_starts_var = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void PiecewiseConstantCoalescent::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    // just update the start times of the intervals
    updateIntervals();
    
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 *
 * @throw RbExpection when no interval start times are specified when the 'SPECIFIED' interval_method is used
 *
 */
void PiecewiseConstantCoalescent::updateIntervals( void )
{
    
    // first check if we want to set the interval ages here!
    if ( interval_method == SPECIFIED )
    {
        if ( interval_starts_var != NULL )
        {
            interval_starts = interval_starts_var->getValue();
        }
        else
        {
            throw RbException("You have to provide the start times of the coalescent skyline intervals if you chose 'SPECIFIED' as the method of choice.");
        }
        
        // clean all the sets
        interval_starts.clear();
        pop_sizes.clear();
        
//        double present_time = process_age->getValue();
        
        std::set<double> event_times;
        
        const std::vector<double>& times = interval_starts_var->getValue();
        for (std::vector<double>::const_iterator it = times.begin(); it != times.end(); ++it)
        {
            event_times.insert( *it );
        }
        
        const std::vector<double> &p = Nes->getValue();
        
        size_t index = p.size()-1;
        
        pop_sizes.push_back( p[index] );
        
        size_t pos = 0;
        for (std::set<double>::iterator it = event_times.begin(); it != event_times.end(); ++it)
        {
            double t = *it;
            
            // add the time to our vector
            interval_starts.push_back( t );
            
            // add the population size at the event t
            pos = size_t( find(times.begin(), times.end(), t) - times.begin() );
            if ( pos != times.size() )
            {
                index = pos;
                pop_sizes.push_back( p[index] );
            }
            else
            {
                pop_sizes.push_back( p[index] );
            }
            
        }
        
        
    }
    else if ( interval_method == EVENTS )
    {
        // we directly use the population size given from the arguments
        pop_sizes = Nes->getValue();
        
        // next, we recompute the starting times of new intervals
        interval_starts = RbVector<double>(Nes->getValue().size()-1, RbConstants::Double::inf);
        
        if ( this->value != NULL )
        {            
            // retrieved the speciation times
            std::vector<double> ages;
            for (size_t i = 0; i < value->getNumberOfInteriorNodes()+1; ++i)
            {
                const TopologyNode& n = value->getInteriorNode( i );
                double a = n.getAge();
                ages.push_back(a);
            }
            // sort the vector of times in ascending order
            std::sort(ages.begin(), ages.end());
            
            size_t num_events_per_interval = size_t( ceil( double(num_taxa-1.0)/Nes->getValue().size()) );
            size_t current_interval = 0;
            size_t current_num_events_in_interval = 0;
            for (size_t i = 0; i < num_taxa-2; ++i)
            {
                ++current_num_events_in_interval;
                if ( current_num_events_in_interval == num_events_per_interval )
                {
                    interval_starts[current_interval] = ages[i];
                    current_num_events_in_interval = 0;
                    ++current_interval;
                }
            }
        }

    }
    
}
