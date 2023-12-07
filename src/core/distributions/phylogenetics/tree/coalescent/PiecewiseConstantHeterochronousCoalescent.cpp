#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>

#include "PiecewiseConstantHeterochronousCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;


/**
 * Default Constructor for the piecewise constant heterochronous coalescent
 *
 * @param N A vector of population sizes for each interval
 * @param i The start time for each interval
 * @param tn A vector of taxon names used during initialization.
 * @param c A vector of clade constraints
 *
 */
PiecewiseConstantHeterochronousCoalescent::PiecewiseConstantHeterochronousCoalescent(const TypedDagNode<RbVector<double> > *N, const TypedDagNode<RbVector<double> > *i, const std::vector<Taxon> &tn, const std::vector<Clade> &c) :
    AbstractCoalescent( tn, c ),
    Nes( N ),
    intervalStarts( i )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( Nes );
    addParameter( intervalStarts );
    
    simulateHeterochronousTree();
}


/** Destructor
 */
PiecewiseConstantHeterochronousCoalescent::~PiecewiseConstantHeterochronousCoalescent()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
PiecewiseConstantHeterochronousCoalescent* PiecewiseConstantHeterochronousCoalescent::clone( void ) const
{
    
    return new PiecewiseConstantHeterochronousCoalescent( *this );
}


/**
 * Compute the log-transformed probability of the current times under the current parameter values.
 *
 * \return    The log-probability density.
 */
double PiecewiseConstantHeterochronousCoalescent::computeLnProbabilityTimes( void ) const
{

    // variable declarations and initialization
    double lnProbTimes = 0;

    const RbVector<double> &popSizes  = Nes->getValue();
    const RbVector<double> &intervals = intervalStarts->getValue();

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
        const TopologyNode& n = value->getTipNode( i );
        double a = n.getAge();
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
        size_t atIntervalStart = 0;
        double nextAge = ages[atAge];
        double nextSerialTime = serialTimes[atSerialTime];
        double nextIntervalStart = intervals[atIntervalStart];

        // create master list of event times and types
        // events are either a sample (lineage size up), coalescence (lineage size down), or theta changepoint (lineage size constant)
        do
        {
            nextAge = ages[atAge];
            if (nextAge <= nextSerialTime && nextAge <= nextIntervalStart) {
                // coalescence
                combinedEventTimes.push_back(nextAge);
                combinedEventTypes.push_back(-1);
                ++atAge;
            } else if (nextSerialTime <= nextAge && nextSerialTime <= nextIntervalStart) {
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
                // theta change
                combinedEventTimes.push_back(nextIntervalStart);
                combinedEventTypes.push_back(0);
                ++atIntervalStart;
                if (atIntervalStart < intervals.size()) {
                    nextIntervalStart = intervals[atIntervalStart];
                } else {
                    nextIntervalStart = RbConstants::Double::inf;
                }
            }
        } while (atAge < ages.size());

    } else {
        size_t atAge = 0;
        size_t atIntervalStart = 0;
        double nextAge = ages[atAge];
        double nextIntervalStart = intervals[atIntervalStart];

        // create master list of event times and types
        // events are either a sample (lineage size up), coalescence (lineage size down), or theta changepoint (lineage size constant)
        do
        {
            nextAge = ages[atAge];
            if (nextIntervalStart <= nextAge) {
                // theta change
                combinedEventTimes.push_back(nextIntervalStart);
                combinedEventTypes.push_back(0);
                ++atIntervalStart;
                if (atIntervalStart < intervals.size()) {
                    nextIntervalStart = intervals[atIntervalStart];
                } else {
                    nextIntervalStart = RbConstants::Double::inf;
                }
            } else {
                // coalescence
                combinedEventTimes.push_back(nextAge);
                combinedEventTypes.push_back(-1);
                ++atAge;
            }
        } while (atAge < ages.size());
    }


    size_t currentInterval = 0;
    size_t j = num_taxaAtPresent;
    double windowStart = 0.0;

    for (size_t i = 0; i < combinedEventTimes.size(); ++i)
    {
        
        double theta = popSizes[currentInterval];
        double nPairs = j * (j-1) / 2.0;
        
        double deltaAge = combinedEventTimes[i] - windowStart;

        if (combinedEventTypes[i] == -1) {
            // coalescence
            lnProbTimes += log( 1.0 / theta ) - nPairs * deltaAge / theta;
            --j;
        } else if (combinedEventTypes[i] == 1) {
            // sampled ancestor
            lnProbTimes -= nPairs * deltaAge / theta ;
            ++j;
        } else {
            // theta change
            lnProbTimes -= nPairs * deltaAge / theta ;
            ++currentInterval;
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
std::vector<double> PiecewiseConstantHeterochronousCoalescent::simulateCoalescentAges( size_t n ) const
{
    // Retreive skyline parameters
    const RbVector<double> &popSizes  = Nes->getValue();
    const RbVector<double> &intervals = intervalStarts->getValue();
    
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
    
    // Put sampling times and pop-size changes into a single vector of event times
    std::vector<double> combinedEventTimes;
    std::vector<double> combinedEventTypes;
    if (num_taxaAtPresent < num_taxa) {
        
        // sort the vector of serial sampling times in ascending order
        std::sort(serialTimes.begin(), serialTimes.end());
        size_t atSerialTime = 0;
        size_t atIntervalStart = 0;
        double nextSerialTime = serialTimes[atSerialTime];
        double nextIntervalStart = intervals[atIntervalStart];
        
        // create master list of event times and types
        // pre-defined events are either a sample (lineage size up) or theta changepoint (lineage size constant)
        // size_t nEvents = 0;
        size_t nTotalEvents = intervals.size() + serialTimes.size();
        for (size_t nEvents = 0; nEvents < nTotalEvents; ++nEvents)
        {
            if (nextIntervalStart <= nextSerialTime) {
                // theta change
                combinedEventTimes.push_back(nextIntervalStart);
                combinedEventTypes.push_back(0.0);
                ++atIntervalStart;
                if (atIntervalStart < intervals.size()) {
                    nextIntervalStart = intervals[atIntervalStart];
                } else {
                    nextIntervalStart = RbConstants::Double::inf;
                }
            } else {
                // serial sample
                combinedEventTimes.push_back(nextSerialTime);
                combinedEventTypes.push_back(1.0);
                ++atSerialTime;
                if (atSerialTime < serialTimes.size()) {
                    nextSerialTime = serialTimes[atSerialTime];
                } else {
                    nextSerialTime = RbConstants::Double::inf;
                }
            }
        }
    } else {
        combinedEventTimes = intervals;
        combinedEventTypes = std::vector<double>(intervals.size(),0.0);
    }
 
    // cap vector with an event at t=infinity
    combinedEventTimes.push_back(RbConstants::Double::inf);
    combinedEventTypes.push_back(RbConstants::Double::inf);

    // allocate the vector for the times
    std::vector<double> coalescentTimes = std::vector<double>(n,0.0);
    
    size_t currentInterval = 0;
    size_t thetaInterval = 0;
    size_t j = num_taxaAtPresent;
    
    // the current age of the process
    double simAge = 0.0;

    // now simulate the ages
    for (size_t i = 0; i < n; ++i)
    {
        
        bool valid = false;
        do
        {
            double nPairs = j * (j-1) / 2.0;
            double theta = popSizes[thetaInterval];
            double lambda = nPairs / theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            simAge += u;
            valid = simAge < combinedEventTimes[currentInterval] && j > 1;
            if ( !valid ) {
                // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                // Alternately, when we cross a serial sampling time or theta window, the number of active lineages changes
                // or the pop size changes, and it is necessary to discard any "excess" time,
                // which is drawn from an incorrect distribution,then we can draw a new time according to
                // the correct number of active lineages.
                // Either we advance or go back, but in both situations we set the time to the current event in combinedEvents.
                simAge = combinedEventTimes[currentInterval];
                if (combinedEventTypes[currentInterval] == 0.0) {
                    // theta change
                    ++thetaInterval;
                } else {
                    // serial sample
                    ++j;
                }
                ++currentInterval;
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
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void PiecewiseConstantHeterochronousCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == Nes)
    {
        Nes = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    else if (oldP == intervalStarts)
    {
        intervalStarts = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
}
