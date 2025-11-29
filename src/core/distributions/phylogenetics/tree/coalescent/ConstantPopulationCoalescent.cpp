#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ConstantPopulationCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Default Constructor
 *
 * @param N a double for the effective population size
 * @param tn a vector of Taxon for the taxon names
 * @param c a vector of clade constraints
 *
 */
ConstantPopulationCoalescent::ConstantPopulationCoalescent(const TypedDagNode<double> *N, const std::vector<Taxon> &tn, const std::vector<Clade> &c) :
    AbstractCoalescent( tn, c ),
    Ne( N )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Ne );
    
    simulateTree();
}


/*
 * Destructor
 */
ConstantPopulationCoalescent::~ConstantPopulationCoalescent()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ConstantPopulationCoalescent* ConstantPopulationCoalescent::clone( void ) const
{
    
    return new ConstantPopulationCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double ConstantPopulationCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double ln_prob = 0;
    double this_ne = Ne->getValue();
    
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
    std::vector<double> serial_tip_ages;
    size_t num_taxa_at_present = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        double a = value->getNode(i).getAge();
        // check if the tip is not contemporaneous
        if ( a > 0.0 )
        {
            serial_tip_ages.push_back(a);
            // remember to decrease the count for the number of taxa at present
            --num_taxa_at_present;
        }
    }
    
    std::vector<double> combined_event_ages;
    std::vector<EVENT_TYPE> combined_event_types;

    // if we have any serially sampled tips
    if (num_taxa_at_present < num_taxa)
    {

        // sort the vector of serial sampling times in ascending order
        std::sort(serial_tip_ages.begin(), serial_tip_ages.end());
        
        size_t at_age = 0;
        size_t at_serial_age = 0;
        double next_age = ages[at_age];
        double next_serial_age = serial_tip_ages[at_serial_age];
        
        // create master list of event times and types
        // events are either a sample (lineage size up), coalescence (lineage size down)
        do
        {
            next_age = ages[at_age];
            if (next_serial_age <= next_age)
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
            else
            {
                // coalescence
                combined_event_ages.push_back( next_age );
                combined_event_types.push_back( COALESCENT );
                ++at_age;
            }
        } while (at_age < ages.size());
        
    }
    else
    {
        combined_event_ages = ages;
        combined_event_types = std::vector<EVENT_TYPE>(ages.size(),COALESCENT);
    }
    
    
    size_t current_num_lineages = num_taxa_at_present;
    // initialize last event age
    double last_event_age = 0.0;
    
    for (size_t i = 0; i < combined_event_ages.size(); ++i)
    {
        double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
        
        double delta_age = combined_event_ages[i] - last_event_age;
        
        if (combined_event_types[i] == SERIAL_SAMPLE)
        {
            // sampled ancestor
            ln_prob -= n_pairs * delta_age / this_ne ;
            ++current_num_lineages;
        }
        else if (combined_event_types[i] == COALESCENT)
        {
            // coalescence probability
            ln_prob += log( 1.0 / this_ne ) - n_pairs * delta_age / this_ne;
            --current_num_lineages;
        }
        else
        {
            throw RbException("Unexpected event type in constant population size coalescent process.");
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
std::vector<double> ConstantPopulationCoalescent::simulateCoalescentAges( size_t n ) const
{
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
    
    size_t at_serial_age = 0;
    if (num_taxa_at_present < num_taxa)
    {
        std::sort(serial_ages.begin(), serial_ages.end());
    }
    
    // now simulate the ages
    
    // allocate the vector for the ages
    std::vector<double> coalescent_ages = std::vector<double>(n,0.0);
    
    // j is the number of active lineages at the current time
    size_t current_num_lineages = num_taxa_at_present;
    double this_ne = Ne->getValue();
    
    // the current age of the process
    double sim_age = 0.0;
    
    // draw a time for each event conditioned on the time of the process
    for (size_t i = 0; i < n; ++i)
    {
        bool is_coalescent_event = false;
        do
        {
            double n_pairs = current_num_lineages * (current_num_lineages-1) / 2.0;
            double lambda = n_pairs / this_ne;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            sim_age += u;
            is_coalescent_event = (at_serial_age >= serial_ages.size() || sim_age < serial_ages[at_serial_age]) && current_num_lineages > 1;
            if ( is_coalescent_event == false )
            {
                // If current_num_lineages is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time, the number of active lineages changes
                // it is necessary to discard any "excess" time, which is drawn from an incorrect distribution
                // then we can draw a new time according to the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current serial sample.
                sim_age = serial_ages[at_serial_age];
                ++at_serial_age;
                ++current_num_lineages;
            }
        } while ( is_coalescent_event == false );
    
        coalescent_ages[i] = sim_age;
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
void ConstantPopulationCoalescent::swapParameterInternal(const DagNode *old_p, const DagNode *new_p)
{
    if (old_p == Ne)
    {
        Ne = static_cast<const TypedDagNode<double>* >( new_p );
    }
}
