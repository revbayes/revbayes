#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>

#include "ConstantPopulationHeterochronousCoalescent.h"
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

ConstantPopulationHeterochronousCoalescent::ConstantPopulationHeterochronousCoalescent(const TypedDagNode<double> *N, const std::vector<Taxon> &tn, const std::vector<Clade> &c) :
    AbstractCoalescent( tn, c ),
    Ne( N )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Ne );
    
    simulateHeterochronousTree();
}



ConstantPopulationHeterochronousCoalescent::~ConstantPopulationHeterochronousCoalescent()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ConstantPopulationHeterochronousCoalescent* ConstantPopulationHeterochronousCoalescent::clone( void ) const
{
    
    return new ConstantPopulationHeterochronousCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double ConstantPopulationHeterochronousCoalescent::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    double theta = Ne->getValue();
    
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
    std::vector<double> serialTimes;
    size_t num_taxaAtPresent = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        double a = value->getNode(i).getAge();
        if ( a > 0.0 ) {
            serialTimes.push_back(a);
            --num_taxaAtPresent;
        }
    }
    
    std::vector<double> combinedEventTimes;
    std::vector<int> combinedEventTypes;

    if (num_taxaAtPresent < num_taxa) {

        // sort the vector of serial sampling times in ascending order
        std::sort(serialTimes.begin(), serialTimes.end());
        
        size_t atAge = 0;
        size_t atSerialTime = 0;
        double nextAge = ages[atAge];
        double nextSerialTime = serialTimes[atSerialTime];
        
        // create master list of event times and types
        // events are either a sample (lineage size up), coalescence (lineage size down), or theta changepoint (lineage size constant)
        do
        {
            nextAge = ages[atAge];
            if (nextSerialTime <= nextAge) {
                // serial sample
                combinedEventTimes.push_back(nextSerialTime);
                combinedEventTypes.push_back(1);
                ++atSerialTime;
                if (atSerialTime < serialTimes.size()) {
                    nextSerialTime = serialTimes[atSerialTime];
                } else {
                    nextSerialTime = RbConstants::Double::inf;
                }
            } else {
                // coalescence
                combinedEventTimes.push_back(nextAge);
                combinedEventTypes.push_back(-1);
                ++atAge;
            }
        } while (atAge < ages.size());
        
    } else {
        combinedEventTimes = ages;
        combinedEventTypes = std::vector<int>(ages.size(),-1);
    }
    
    
    size_t j = num_taxaAtPresent;
    double windowStart = 0.0;
    
    for (size_t i = 0; i < combinedEventTimes.size(); ++i)
    {
        double nPairs = j * (j-1) / 2.0;
        
        double deltaAge = combinedEventTimes[i] - windowStart;
        
        if (combinedEventTypes[i] == 1) {
            // sampled ancestor
            lnProbTimes -= nPairs * deltaAge / theta ;
            ++j;
        } else {
            // coalescence
            lnProbTimes += log( 1.0 / theta ) - nPairs * deltaAge / theta;
            --j;
        }
        
        windowStart = combinedEventTimes[i];
    }
    
    return lnProbTimes;
    
}

/**
 * Simulate new coalescent times.
 *
 * \param[in]    n      The number of coalescent events to simulate.
 *
 * \return    A vector of the simulated coalescent times.
 */
std::vector<double> ConstantPopulationHeterochronousCoalescent::simulateCoalescentAges( size_t n ) const
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serialTimes;
    size_t num_taxaAtPresent = 0;
    for (size_t i = 0; i < num_taxa; ++i)
    {
        double a = taxa[i].getAge();
        if ( a > 0.0 ) {
            serialTimes.push_back(a);
        } else {
            ++num_taxaAtPresent;
        }
    }
    
    size_t atSerialTime = 0;
    if (num_taxaAtPresent == num_taxa) {
        serialTimes.push_back(RbConstants::Double::inf);
    } else {
        std::sort(serialTimes.begin(), serialTimes.end());
        serialTimes.push_back(RbConstants::Double::inf);
    }
    
    // now simulate the ages
    
    // allocate the vector for the times
    std::vector<double> coalescentTimes = std::vector<double>(n,0.0);
    
    // j is the number of active lineages at the current time
    size_t j = num_taxaAtPresent;
    double theta = Ne->getValue();
    
    // the current age of the process
    double simAge = 0.0;
    
    // draw a time for each speciation event condition on the time of the process
    for (size_t i = 0; i < n; ++i)
    {
        bool valid = false;
        do
        {
            double nPairs = j * (j-1) / 2.0;
            double lambda = nPairs / theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            simAge += u;
            valid = simAge < serialTimes[atSerialTime] && j > 1;
            if ( !valid ) {
                // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time, the number of active lineages changes
                // it is necessary to discard any "excess" time, which is drawn from an incorrect distribution
                // then we can draw a new time according to the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current serial sample.
                simAge = serialTimes[atSerialTime];
                ++atSerialTime;
                ++j;
            }
        } while ( !valid );
    
        coalescentTimes[i] = simAge;
        --j;
        
    }
        
    return coalescentTimes;
}

/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void ConstantPopulationHeterochronousCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == Ne)
    {
        Ne = static_cast<const TypedDagNode<double>* >( newP );
    }
}
