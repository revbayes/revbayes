#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
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
 * \param    t         The starting tree if we want to avoid simulating trees. 
 */
AbstractBirthDeathProcess::AbstractBirthDeathProcess(const TypedDagNode<double> *ra, const std::string &cdt, const std::vector<Taxon> &tn, bool uo, Tree* t )
    : AbstractRootedTreeDistribution( ra, tn, uo, t ),
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
        size_t n_taxa_present = getNumberOfTaxaAtPresent();
        
        lnProbTimes = -lnProbNumTaxa( n_taxa_present, 0, present_time, true );
    }
    
    // multiply the probability of a descendant of the initial species
    lnProbTimes += computeLnProbabilityTimes();
    
    return lnProbTimes;
}


size_t AbstractBirthDeathProcess::getNumberOfTaxaAtPresent( void ) const
{
    
    size_t num_taxa_present = value->getNumberOfTips();
    
    return num_taxa_present;
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
 * @param alwaysReturn whether the simulation can return times which are not valid draws from the distribution (for initial values)
 * @return vector of simulated divergence times
 **/
std::vector<double> AbstractBirthDeathProcess::simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const
{

    std::vector<double> times(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        double t = 0.0;
        int ntries = 0;
        
        do { 
            t = simulateDivergenceTime(origin, present);
            ntries ++;
        } while (t < min && ntries < 1000);

        if(t < min) {
            if(!alwaysReturn) throw RbException("No valid clade age found for simulated tree after 1000 tries.");

            auto rng = GLOBAL_RNG;
            double u = rng->uniform01();
            t = min + u*(origin-min);
        }

        times[i] = t;
    }
    
    // finally sort the times
    std::sort(times.begin(), times.end());
    
    return times;
}

