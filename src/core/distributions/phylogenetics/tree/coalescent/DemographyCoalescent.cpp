#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "DemographyCoalescent.h"
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
DemographyCoalescent::DemographyCoalescent(const TypedDagNode< RbVector<double> > *iv, const RbVector< DemographicFunction > &df, const std::vector<Taxon> &tn, const std::vector<Clade> &c) : AbstractCoalescent( tn, c ),
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

DemographyCoalescent::~DemographyCoalescent()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
DemographyCoalescent* DemographyCoalescent::clone( void ) const
{
    
    return new DemographyCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double DemographyCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double ln_prob = 0;
    
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
    std::vector<double> serial_ages;
    size_t num_taxa_at_present = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        double a = value->getTipNode(i).getAge();
        if ( a > 0.0 )
        {
            serial_ages.push_back(a);
            --num_taxa_at_present;
        }
    }
    
    // retrieve change times
    const RbVector<double> &change_ages = intervals->getValue();
    
    std::vector<double> combined_event_ages;
    std::vector<EVENT_TYPE> combined_event_types;
    
    bool heterochronous = num_taxa_at_present < num_taxa;
    
    if ( heterochronous == true )
    {
        // sort the vector of serial sampling times in ascending order
        std::sort(serial_ages.begin(), serial_ages.end());
    }

    size_t index_age = 0;
    size_t index_serial_age = 0;
    size_t index_demographic_function_change_point = 0;
    double next_age = ages[index_age];
    double next_serial_age = RbConstants::Double::nan;
    if ( heterochronous == true )
    {
        next_serial_age = serial_ages[index_serial_age];
    }
    double next_df_change_age = RbConstants::Double::inf;
    if ( index_demographic_function_change_point < change_ages.size() )
    {
        next_df_change_age = change_ages[index_demographic_function_change_point];
    }
    // create master list of event times and types
    // events are either a sample (lineage size up), coalescence (lineage size down), or Ne changepoint (lineage size constant)
    do
    {
        next_age = ages[index_age];
        if ( heterochronous == true && next_serial_age <= next_age && next_serial_age <= next_df_change_age )
        {
            // serial sample
            combined_event_ages.push_back( next_serial_age );
            combined_event_types.push_back( SERIAL_SAMPLE );
            ++index_serial_age;
            if (index_serial_age < serial_ages.size())
            {
                next_serial_age = serial_ages[index_serial_age];
            }
            else
            {
                next_serial_age = RbConstants::Double::inf;
            }
        }
        else if ( next_df_change_age <= next_age )
        {
            // change of demographic function
            combined_event_ages.push_back(next_df_change_age);
            combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
            ++index_demographic_function_change_point;
            if ( index_demographic_function_change_point < change_ages.size() )
            {
                next_df_change_age = change_ages[index_demographic_function_change_point];
            }
            else
            {
                next_df_change_age = RbConstants::Double::inf;
            }
        }
        else
        {
            // coalescence
            combined_event_ages.push_back(next_age);
            combined_event_types.push_back(COALESCENT);
            ++index_age;
        }
    } while (index_age < ages.size());
        
    
    
    size_t current_num_lineages = num_taxa_at_present;
    double last_event_age = 0.0;
    size_t index_demographic_function = 0;
    const DemographicFunction *current_demographic_function = &demographies[index_demographic_function];
    
    
    for (size_t i = 0; i < combined_event_ages.size(); ++i)
    {
        double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
        double interval_area = current_demographic_function->getIntegral(last_event_age, combined_event_ages[i]);
        
        // add log probability that nothing happens until the next event
        ln_prob -= n_pairs * interval_area;
        
        // handle events accordingly
        if (combined_event_types[i] == SERIAL_SAMPLE )
        {
            // sampled ancestor
            ++current_num_lineages;
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
            double ne_at_coal_age = current_demographic_function->getDemographic(combined_event_ages[i]);
            ln_prob -= log( ne_at_coal_age );
            --current_num_lineages;
        }
        
        last_event_age = combined_event_ages[i];
    }
    
    return ln_prob;
    
}

/**
 * Simulate new coalescent times.
 *
 * \param[in]    n      The number of coalescent events to simulate.
 *
 * \return    A vector of the simulated coalescent times.
 */
std::vector<double> DemographyCoalescent::simulateCoalescentAges( size_t n ) const
{
    const RbVector<double> &change_ages = intervals->getValue();

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serial_ages;
    size_t num_taxa_at_present = 0;
    for (size_t i = 0; i < num_taxa; ++i)
    {
        double a = taxa[i].getAge();
        if ( a > 0.0 )
        {
            serial_ages.push_back(a);
        }
        else
        {
            ++num_taxa_at_present;
        }
    }

    // Put sampling times and demographic function changes into a single vector of event times
    std::vector<double> combined_event_ages;
    std::vector<double> combined_event_types;
    if (num_taxa_at_present < num_taxa)
    {
        // sort the vector of serial sampling times in ascending order
        std::sort(serial_ages.begin(), serial_ages.end());
        size_t atSerialTime = 0;
        size_t atIntervalStart = 0;
        double nextSerialTime = serial_ages[atSerialTime];
        double nextIntervalStart = change_ages[atIntervalStart];
        
        // create master list of event times and types
        // pre-defined events are either a sample (lineage size up) or demographic function changepoint (lineage size constant)
        size_t nTotalEvents = change_ages.size() + serial_ages.size();
        for (size_t nEvents = 0; nEvents < nTotalEvents; ++nEvents)
        {
            if (nextIntervalStart <= nextSerialTime)
            {
                // demographic function change
                combined_event_ages.push_back(nextIntervalStart);
                combined_event_types.push_back( DEMOGRAPHIC_MODEL_CHANGE );
                ++atIntervalStart;
                if (atIntervalStart < change_ages.size())
                {
                    nextIntervalStart = change_ages[atIntervalStart];
                }
                else
                {
                    nextIntervalStart = RbConstants::Double::inf;
                }
            }
            else
            {
                // serial sample
                combined_event_ages.push_back(nextSerialTime);
                combined_event_types.push_back( SERIAL_SAMPLE );
                ++atSerialTime;
                if (atSerialTime < serial_ages.size())
                {
                    nextSerialTime = serial_ages[atSerialTime];
                }
                else
                {
                    nextSerialTime = RbConstants::Double::inf;
                }
            }
        }
    }
    else
    {
        combined_event_ages = change_ages;
        combined_event_types = std::vector<double>(change_ages.size(),DEMOGRAPHIC_MODEL_CHANGE);
    }
 
    // cap vector with an event at t=infinity
    combined_event_ages.push_back(RbConstants::Double::inf);
    combined_event_types.push_back(RbConstants::Double::inf);

    
    // now simulate the ages
    
    // allocate the vector for the times
    std::vector<double> coalescent_ages = std::vector<double>(n,0.0);
    
    // j is the number of active lineages at the current time
    size_t j = num_taxa_at_present;
    size_t currentInterval = 0;
    // size_t index_serial_age = 0;
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

            valid = (sim_age < combined_event_ages[currentInterval] && j > 1) && waitingTime > 0;
            if ( valid == false )
            {
                // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time or Ne window, the number of active lineages changes
                // or the pop size changes, and it is necessary to discard any "excess" time,
                // which is drawn from an incorrect distribution,then we can draw a new time according to
                // the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current event in combinedEvents.
                sim_age = combined_event_ages[currentInterval];
                if (combined_event_types[currentInterval] == DEMOGRAPHIC_MODEL_CHANGE)
                {
                    // demographic function change
                    ++index_demographic_function;
                    current_demographic_function = &demographies[index_demographic_function];
                }
                else
                {
                    // serial sample
                    ++j;
                }
                ++currentInterval;
            }
            
        } while ( valid == false );
        
        coalescent_ages[i] = sim_age;
        --j;
        
    }
    
    return coalescent_ages;
}

/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void DemographyCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
