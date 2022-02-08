#include <stddef.h>
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
    interval_change_points_var( i ),
    interval_method( meth )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Nes );
    addParameter( interval_change_points_var );
    
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
    
    // initialize the current interval index at 0
    size_t current_interval = 0;
    
    size_t current_num_lineages = num_taxa;
    double last_event_age = 0.0;
    for (size_t i = 0; i < coalescent_ages.size(); ++i)
    {
        double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;

        // initialize theta
        double theta = 0.0;
        
        // initialize the waiting time variable
        double delta_age = coalescent_ages[i];
        bool was_coalescent_event = false;
        do
        {
            theta = pop_sizes[current_interval];
            double next_coal_age = coalescent_ages[i];
            double next_event_age = next_coal_age;
            if ( current_interval < interval_change_points.size() )
            {
                next_event_age = (next_coal_age > interval_change_points[current_interval]) ? interval_change_points[current_interval] : next_coal_age;
            }
            
            delta_age = next_event_age - last_event_age;
            was_coalescent_event = current_interval >= interval_change_points.size() || next_coal_age <= interval_change_points[current_interval];
            if ( was_coalescent_event == false )
            {
                ++current_interval;
                last_event_age = next_event_age;
                
                // we multiply with the waiting time that no coalescent event happened until the next interval change point
                ln_prob -= n_pairs * delta_age / theta ;
            }
            
        } while ( was_coalescent_event == false );
        
        // we multiply with the probability of a coalescent event
        ln_prob += log( 1.0 / theta ) - n_pairs * delta_age / theta;
        
        // decrement the current number of lineage counter
        --current_num_lineages;
    }
    
    return ln_prob;
}


void PiecewiseConstantCoalescent::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "getIntervalAges" )
    {
        
        rv = interval_change_points;
        
    }
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

    // Sebastian: This is currently redudant because we update the intervals each time when we compute the probability
    // just re-update the start times of the intervals
//    updateIntervals();
    
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
        updateIntervals();
    }
    else
    {
        interval_change_points = RbVector<double>(Nes->getValue().size()-1, RbConstants::Double::inf);
    }
    size_t num_events_per_interval = size_t( floor( double(n)/Nes->getValue().size()) );
    size_t current_num_events_in_interval = 0;
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // now simulate the ages
    
    // allocate the vector for the times
    std::vector<double> coalescent_ages = std::vector<double>(n,0.0);
    
    const RbVector<double> &pop_sizes  = Nes->getValue();
    size_t current_interval = 0;
    size_t current_num_lineages = num_taxa;
    // draw a time for each speciation event condition on the time of the process
    for (size_t i = 0; i < n; ++i)
    {
        double prev_coalescent_age = 0.0;
        if ( i > 0 )
        {
            prev_coalescent_age = coalescent_ages[i-1];
        }
        
        double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
        
        double sim_age = 0.0;
        bool was_coalescent_event = false;
        do
        {
            double theta = 1.0 / (pop_sizes[current_interval]);
            double lambda = n_pairs * theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            sim_age = prev_coalescent_age + u;
            was_coalescent_event = current_interval >= interval_change_points.size() || sim_age < interval_change_points[current_interval];
            if ( was_coalescent_event == false )
            {
                prev_coalescent_age = interval_change_points[current_interval];
                ++current_interval;
            }
            
        } while ( was_coalescent_event == false );
        
        coalescent_ages[i] = sim_age;
        
        ++current_num_events_in_interval;
        if ( interval_method == EVENTS && current_num_events_in_interval == num_events_per_interval && current_interval < interval_change_points.size() )
        {
            interval_change_points[current_interval] = sim_age;
            current_num_events_in_interval = 0;
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
void PiecewiseConstantCoalescent::swapParameterInternal(const DagNode *old_p, const DagNode *new_p)
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
void PiecewiseConstantCoalescent::touchSpecialization(const DagNode *affecter, bool touchAll)
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
void PiecewiseConstantCoalescent::updateIntervals( void ) const
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
            
            size_t num_events_per_interval = size_t( floor( double(num_taxa-1.0)/Nes->getValue().size()) );
            size_t current_interval = 0;
            size_t current_num_events_in_interval = 0;
            for (size_t i = 0; i < num_taxa-2; ++i)
            {
                ++current_num_events_in_interval;
                if ( current_num_events_in_interval == num_events_per_interval )
                {
                    interval_change_points[current_interval] = coalescent_ages[i];
                    current_num_events_in_interval = 0;
                    ++current_interval;
                }
            }
        }

    }
    
}
