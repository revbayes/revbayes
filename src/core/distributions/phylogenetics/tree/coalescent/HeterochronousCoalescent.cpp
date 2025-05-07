#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "HeterochronousCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "DagNode.h"
#include "DemographicFunction.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;


/**
 * Default Constructor for Heterochronus Coalescent
 *
 * @param[in] iv The start times of intervals where demographic dynamics can differ
 * @param[in] df A vector specifying the demographic dynamics for a given interval.
 * @param[in] tn        A vector of taxon names used during initialization.
 * @param[in] c         A vector of clade constraints.
 */
HeterochronousCoalescent::HeterochronousCoalescent(const TypedDagNode< RbVector<double> > *iv, const RbVector< DemographicFunction > &df, const std::vector<Taxon> &tn, const std::vector<Clade> &c) : AbstractCoalescent( tn, c ),
    intervals( iv ),
    demographies( df )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( intervals );
    
    for (size_t i=0; i<demographies.size(); ++i)
    {
        const std::vector<const DagNode*> &pars = demographies[i].getDagNodes();
        for (size_t j=0; j<pars.size(); ++j)
        {
            addParameter( pars[j] );
        }
    }
    
    simulateHeterochronousTree();
}

HeterochronousCoalescent::~HeterochronousCoalescent()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
HeterochronousCoalescent* HeterochronousCoalescent::clone( void ) const
{
    
    return new HeterochronousCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double HeterochronousCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double ln_prob_times = 0;
    
    // retrieve the coalescence times
    std::vector<double> ages;
    for (size_t i = 0; i < value->getNumberOfInteriorNodes()+1; ++i)
    {
        const TopologyNode& n = value->getInteriorNode( i );
        double a = n.getAge();
        ages.push_back(a);
    }
    // sort the vector of coalescence times in ascending order
    std::sort(ages.begin(), ages.end());
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serial_times;
    size_t num_taxa_at_present = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        double a = value->getTipNode(i).getAge();
        if ( a > 0.0 )
        {
            serial_times.push_back(a);
            --num_taxa_at_present;
        }
    }
    
    // retrieve change times
    const RbVector<double> &change_times = intervals->getValue();
    
    std::vector<double> combined_event_times;
    std::vector<EVENT_TYPE> combined_event_types;
    
    bool heterochronous = num_taxa_at_present < num_taxa;
    
    if ( heterochronous == true )
    {
        // sort the vector of serial sampling times in ascending order
        std::sort(serial_times.begin(), serial_times.end());
    }

    size_t index_age = 0;
    size_t index_serial_time = 0;
    size_t index_demographic_function_change_point = 0;
    double next_age = ages[index_age];
    double next_serial_time = RbConstants::Double::nan;
    if ( heterochronous == true )
    {
        next_serial_time = serial_times[index_serial_time];
    }
    double next_df_change_time = RbConstants::Double::inf;
    if ( index_demographic_function_change_point < change_times.size() )
    {
        next_df_change_time = change_times[index_demographic_function_change_point];
    }
    // create master list of event times and types
    // events are either a sample (lineage size up), coalescence (lineage size down), or theta changepoint (lineage size constant)
    do
    {
        next_age = ages[index_age];
        if ( heterochronous == true && next_serial_time <= next_age && next_serial_time <= next_df_change_time )
        {
            // serial sample
            combined_event_times.push_back(next_serial_time);
            combined_event_types.push_back( SERIAL_SAMPLE );
            ++index_serial_time;
            if (index_serial_time < serial_times.size())
            {
                next_serial_time = serial_times[index_serial_time];
            }
            else
            {
                next_serial_time = RbConstants::Double::inf;
            }
        }
        else if ( next_df_change_time <= next_age )
        {
            // change of demographic function
            combined_event_times.push_back(next_df_change_time);
            combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
            ++index_demographic_function_change_point;
            if ( index_demographic_function_change_point < change_times.size() )
            {
                next_df_change_time = change_times[index_demographic_function_change_point];
            }
            else
            {
                next_df_change_time = RbConstants::Double::inf;
            }
        }
        else
        {
            // coalescence
            combined_event_times.push_back(next_age);
            combined_event_types.push_back(COALESCENT);
            ++index_age;
        }
    } while (index_age < ages.size());
        
    
    
    size_t j = num_taxa_at_present;
    double window_start = 0.0;
    size_t index_demographic_function = 0;
    const DemographicFunction *current_demographic_function = &demographies[index_demographic_function];
    
    
    for (size_t i = 0; i < combined_event_times.size(); ++i)
    {
        double n_pairs = j * (j-1) / 2.0;
        double interval_area = current_demographic_function->getIntegral(window_start, combined_event_times[i]);
        
        // add log probability that nothing happens until the next event
        ln_prob_times -= n_pairs * interval_area;
        
        // handle events accordingly
        if (combined_event_types[i] == SERIAL_SAMPLE )
        {
            // sampled ancestor
            ++j;
        }
        else if ( combined_event_types[i] == DEMOGRAPHIC_MODEL_CHANGE )
        {
            // change of the demographic function
            ++index_demographic_function;
            if ( index_demographic_function > demographies.size() )
            {
                throw RbException("Problem occurred in coalescent process with demographic functions: We tried to access a demographic function outside the vector.");
            }
            current_demographic_function = &demographies[index_demographic_function];
        }
        else
        {
            // coalescence
            double theta_at_coal_time = current_demographic_function->getDemographic(combined_event_times[i]);
            ln_prob_times -= log( theta_at_coal_time );
            --j;
        }
        
        window_start = combined_event_times[i];
    }
    
    return ln_prob_times;
    
}

/**
 * Simulate new coalescent times.
 *
 * \param[in]    n      The number of coalescent events to simulate.
 *
 * \return    A vector of the simulated coalescent times.
 */
std::vector<double> HeterochronousCoalescent::simulateCoalescentAges( size_t n ) const
{
    const RbVector<double> &change_times = intervals->getValue();

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serial_times;
    size_t num_taxa_at_present = 0;
    for (size_t i = 0; i < num_taxa; ++i)
    {
        double a = taxa[i].getAge();
        if ( a > 0.0 )
        {
            serial_times.push_back(a);
        }
        else
        {
            ++num_taxa_at_present;
        }
    }

    // Put sampling times and demographic function changes into a single vector of event times
    std::vector<double> combinedEventTimes;
    std::vector<double> combinedEventTypes;
    if (num_taxa_at_present < num_taxa) {
        // sort the vector of serial sampling times in ascending order
        std::sort(serial_times.begin(), serial_times.end());
        size_t atSerialTime = 0;
        size_t atIntervalStart = 0;
        double nextSerialTime = serial_times[atSerialTime];
        double nextIntervalStart = change_times[atIntervalStart];
        
        // create master list of event times and types
        // pre-defined events are either a sample (lineage size up) or demographic function changepoint (lineage size constant)
        size_t nTotalEvents = change_times.size() + serial_times.size();
        for (size_t nEvents = 0; nEvents < nTotalEvents; ++nEvents)
        {
            if (nextIntervalStart <= nextSerialTime) {
                // demographic function change
                combinedEventTimes.push_back(nextIntervalStart);
                combinedEventTypes.push_back(0.0);
                ++atIntervalStart;
                if (atIntervalStart < change_times.size()) {
                    nextIntervalStart = change_times[atIntervalStart];
                } else {
                    nextIntervalStart = RbConstants::Double::inf;
                }
            } else {
                // serial sample
                combinedEventTimes.push_back(nextSerialTime);
                combinedEventTypes.push_back(1.0);
                ++atSerialTime;
                if (atSerialTime < serial_times.size()) {
                    nextSerialTime = serial_times[atSerialTime];
                } else {
                    nextSerialTime = RbConstants::Double::inf;
                }
            }
        }
    } else {
        combinedEventTimes = change_times;
        combinedEventTypes = std::vector<double>(change_times.size(),0.0);
    }
 
    // cap vector with an event at t=infinity
    combinedEventTimes.push_back(RbConstants::Double::inf);
    combinedEventTypes.push_back(RbConstants::Double::inf);

    
    // now simulate the ages
    
    // allocate the vector for the times
    std::vector<double> coalescent_times = std::vector<double>(n,0.0);
    
    // j is the number of active lineages at the current time
    size_t j = num_taxa_at_present;
    size_t currentInterval = 0;
    // size_t index_serial_time = 0;
    size_t index_demographic_function = 0;

    const DemographicFunction *current_demographic_function = &demographies[index_demographic_function];

    // the current age of the process
    double sim_age = 0.0;
    
    // draw a time for each speciation event condition on the time of the process
    for (size_t i = 0; i < n; ++i)
    {
        bool valid = false;
        do
        {
            double nPairs = j * (j-1) / 2.0;
            double lambda = RbStatistics::Exponential::rv( nPairs, *rng);
            double waitingTime = current_demographic_function->getWaitingTime(sim_age, lambda);
            sim_age += waitingTime;

            valid = (sim_age < combinedEventTimes[currentInterval] && j > 1) && waitingTime > 0;
            if ( valid == false ) {
                // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time or theta window, the number of active lineages changes
                // or the pop size changes, and it is necessary to discard any "excess" time,
                // which is drawn from an incorrect distribution,then we can draw a new time according to
                // the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current event in combinedEvents.
                sim_age = combinedEventTimes[currentInterval];
                if (combinedEventTypes[currentInterval] == 0.0) {
                    // demographic function change
                    ++index_demographic_function;
                    current_demographic_function = &demographies[index_demographic_function];
                } else {
                    // serial sample
                    ++j;
                }
                ++currentInterval;
            }
            
        } while ( valid == false );
        
        coalescent_times[i] = sim_age;
        --j;
        
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
void HeterochronousCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    bool found = false;
    if ( oldP == intervals )
    {
        intervals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
        found = true;
    }
    else
    {
        for (size_t i=0; i<demographies.size(); ++i)
        {
            
            try {
                demographies[i].swapNode(oldP, newP);
                // if the statement succeeded and didn't throw an error, then the distribution had this parameter
                found = true;
            }
            catch (RbException &e)
            {
                // do nothing because we actually do not know who had the parameter
            }
        
        }
    }
    
    
    if ( found == false )
    {
        throw RbException("Could not find the distribution parameter to be swapped for the demographic function coalescent process: " + oldP->getName() + " to " + newP->getName()) ;
    }
    
}
