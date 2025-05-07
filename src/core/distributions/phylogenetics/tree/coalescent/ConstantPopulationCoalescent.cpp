#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ConstantPopulationCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "AbstractCoalescent.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;


/**
 * Default Constructor
 *
 * @param N a double for the effective population size
 * @param tn a vector of Taxon for the taxon names
 * @param c a vector of clade constraints
 *
 */
ConstantPopulationCoalescent::ConstantPopulationCoalescent(const TypedDagNode<double> *N, const std::vector<Taxon> &tn, const std::vector<Clade> &c) : AbstractCoalescent( tn, c ),
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
    double ln_prob_times = 0;
    
    
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
    
    for (size_t i = 0; i < ages.size(); ++i) 
    {
//        size_t j = num_taxa - i;
//        double theta = 1.0 / (2.0*Ne->getValue());
//        double nPairs = j * (j-1) / 2.0;
//        lnProbTimes += log( nPairs * 2.0 / theta ) - nPairs * 2.0 / theta * ages[i] ;
        
        
        size_t j = num_taxa - i;
        double theta = Ne->getValue();
        double nPairs = j * (j-1) / 2.0;
        
        double prevCoalescentTime = 0.0;
        if ( i > 0 )
        {
            prevCoalescentTime = ages[i-1];
        }
        
        double deltaAge = ages[i] - prevCoalescentTime;

        // CoalescentMCMC
//        lnProbTimes += log( nPairs / theta ) - 2 * nPairs * deltaAge / theta;
        // BEAST:
//        lnProbTimes += log( 1.0 / theta ) - nPairs * deltaAge / theta ;
        // RevBayes"
//        lnProbTimes += log( nPairs / theta ) - nPairs * deltaAge / theta;
        
        
        ln_prob_times += log( 1.0 / theta ) - nPairs * deltaAge / theta;
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
std::vector<double> ConstantPopulationCoalescent::simulateCoalescentAges( size_t n ) const
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // allocate the vector for the times
    std::vector<double> coalescent_times = std::vector<double>(n,0.0);
    
//    double theta = 1.0 / (2.0*Ne->getValue());
    double theta = 1.0 / (Ne->getValue());
    // draw a time for each speciation event condition on the time of the process
	for (size_t i = 0; i < n; ++i)
    {
        double prev_coalescent_time = 0.0;
        if ( i > 0 ) 
        {
            prev_coalescent_time = coalescent_times[i-1];
        }
        
        size_t j = num_taxa - i;
        double nPairs = j * (j-1) / 2.0;
        double lambda = nPairs * theta;
        double u = RbStatistics::Exponential::rv( lambda, *rng);
		coalescent_times[i] = prev_coalescent_time + u;
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
void ConstantPopulationCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == Ne)
    {
        Ne = static_cast<const TypedDagNode<double>* >( newP );
    }
}
