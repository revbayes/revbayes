#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
#include "AbstractRootedTreeDistribution.h"
#include "Tree.h"

namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 *
 * \param    ra        Origin or root age of the process.
 * \param    cdt       The condition of the process (time/survival/nTaxa)
 * \param    tn        Taxon names used during initialization.
 * \param    uo        If true ra is the origin time otherwise the root age of the process.
 */
AbstractBirthDeathProcess::AbstractBirthDeathProcess(const TypedDagNode<double> *ra, const std::string &cdt, const std::vector<Taxon> &tn, bool uo )
    : AbstractRootedTreeDistribution( ra, tn, uo ),
    condition( cdt )
{}


AbstractBirthDeathProcess::~AbstractBirthDeathProcess(void)
{}


/**
 * Prepare and run the computation of the full log-likelihood of the divergence times.
 *
 * @return log-probability
 */
double AbstractBirthDeathProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    // prepare the probability computation
    prepareProbComputation();
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    
    // present time
    double present_time = value->getRoot().getAge();
    
    // what do we condition on?
    // did we condition on survival?
    if ( condition == "survival" )
    {
        lnProbTimes = - log( pSurvival(0, present_time) );
        
        // we started at the root thus we square the survival prob
        lnProbTimes *= 2.0;
        
    }
    else if ( condition == "nTaxa" )
    {
        lnProbTimes = -lnProbNumTaxa( value->getNumberOfTips(), 0, present_time, true );
    }
    
    // multiply the probability of a descendant of the initial species
    lnProbTimes += computeLnProbabilityTimes();
    
    return lnProbTimes;
}


/** The birth death density is derived for a (ranked) unlabeled oriented tree
 * so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n-1} / n!
 * see Gernhard 2008
 *
 * @return log of the computed factor
 **/
double AbstractBirthDeathProcess::lnProbTreeShape(void) const
{
    return (value->getNumberOfTips() - 1) * RbConstants::LN2 - RbMath::lnFactorial( (int)value->getNumberOfTips() );
}


/**
 * Prepare the probability computation. Here we can pre-calculate some values for more
 * efficient probability calculation. The derived classes may want to do something.
 */
void AbstractBirthDeathProcess::prepareProbComputation( void ) const
{}


/**
 * Simulate divergence times under a birth-death process.
 *
 * @param n number of times to simulate
 * @param origin start time of the process
 * @param present stop time of the process
 * @param min minimum value of the simulated times
 * @return vector of simulated divergence times
 **/
std::vector<double> AbstractBirthDeathProcess::simulateDivergenceTimes(size_t n, double origin, double present, double min) const
{

    std::vector<double> times(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        double t = 0.0;

        do {
            t = simulateDivergenceTime(origin, present);
        } while ( t < min );

        times[i] = t;
    }
    
    // finally sort the times
    std::sort(times.begin(), times.end());
    
    return times;
}


TopologyNode* AbstractBirthDeathProcess::simulateOneHiddenClade(double t) {
     throw RbException("AbstractBirthDeathProcess::simulateOneHiddenClade() not implemented.");
}

/**
 * Simulate hidden clades upon the reconstructed tree under a birth-death process.
 **/

void AbstractBirthDeathProcess::simulateHiddenClades(void) {
    
    throw RbException("AbstractBirthDeathProcess::simulateHiddenClades() not implemented.");

}
