#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iosfwd>
#include <vector>

#include "DistributionExponential.h"
#include "DiversityDependentPureBirthProcess.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "AbstractBirthDeathProcess.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;


/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    o      Time of the origin/present/length of the process.
 * \param[in]    s      Initial speciation rate (lambda_0).
 * \param[in]    k      Carrying capacity.
 * \param[in]    cdt    Condition of the process (none/survival/#Taxa).
 * \param[in]    tn     Taxa.
 * \param[in]    c      Clades conditioned to be present.
 */
DiversityDependentPureBirthProcess::DiversityDependentPureBirthProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<std::int64_t> *k,
                                                                       const std::string &cdt, const std::vector<Taxon> &tn) : AbstractBirthDeathProcess( ra, cdt, tn, false, NULL ),
        initialSpeciation( s ), 
        capacity( k ) 
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( initialSpeciation );
    addParameter( capacity );
    
    simulateTree(true);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
DiversityDependentPureBirthProcess* DiversityDependentPureBirthProcess::clone( void ) const
{
    
    return new DiversityDependentPureBirthProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double DiversityDependentPureBirthProcess::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double lnProbTimes = 0;
    
    // present time 
    double tipAge = value->getTipNode(0).getAge();
    
    // present time
    double ra = value->getRoot().getAge();
    double presentTime = ra;
    
    // test that the time of the process is larger or equal to the present time
    if ( tipAge > presentTime )
    {
        return RbConstants::Double::neginf;
    }
    
    // retrieved the speciation times
    recomputeDivergenceTimesSinceOrigin();
    
    int n = 1;
    double b = initialSpeciation->getValue();
    int k = (int)capacity->getValue();
    double lastTime = 0.0;
    double speciationRate, timeInterval;
    for (size_t i = 1; i < value->getNumberOfTips()-1; ++i)
    {
        if ( RbMath::isNan(lnProbTimes) ||
            lnProbTimes == RbConstants::Double::inf || 
            lnProbTimes == RbConstants::Double::neginf ) 
        {
            return RbConstants::Double::nan;
        }
        
        speciationRate = (1.0 - double(n)/k) * b ;
        timeInterval = divergence_times[i] - lastTime;
        lastTime = divergence_times[i];
        
        lnProbTimes += log(speciationRate) - double(n) * speciationRate * timeInterval;
        ++n;
    }
    speciationRate = double(n) * (1.0 - double(n)/k) * b ;
    timeInterval = presentTime - lastTime;
    lnProbTimes -= speciationRate * timeInterval;
    
    return lnProbTimes;
    
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 * The general equation for the homogeneous birth-death process is
 *    P(N(T)>0|N(t)=1) = (1+\int{mu*exp[rateIntegral(t,s)]ds})^{-1}
 * Because this process has no extinction, the survival probability is simplified to
 *    P(N(T)>0|N(t)=1) = 1
 * 
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Speciation rate at time t. 
 */
double DiversityDependentPureBirthProcess::pSurvival(double start, double end) const
{
    
    return 1.0;
}


/**
 * Simulate new speciation times.
 */
double DiversityDependentPureBirthProcess::simulateDivergenceTime(double origin, double present) const
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the parameters
    double lambda = initialSpeciation->getValue();
    double k = capacity->getValue();
    int n = 1;
    
    // \todo
    // draw the final event
    // this is not until actually an event happened but a uniform time before the next species would have been sampled.
    
    double rate = fmax( 1.0 - ((n+3)/k), 1E-8 ) * lambda;
    double t = origin;
    
    do {
        t = RbStatistics::Exponential::rv(rate, *rng);
    } while ( present + t > origin );
    
    
    
    return present + t;
    
//    double lastEvent = t * rng->uniform01();
//    
//    std::vector<double> *times = new std::vector<double>(n,0.0);
//    (*times)[n-1] = lastEvent;
//    for (size_t i = 1; i < n; i++ )
//    {
//        rate = ( 1.0 - ((n-i+2)/k) ) * lambda;
//        t = lastEvent + RbStatistics::Exponential::rv(rate, *rng);
//        lastEvent = t;
//        (*times)[n-i-1] = t;
//    }
//    
//    rate = ( 1.0 - (2/k) ) * lambda;
//    lastEvent += RbStatistics::Exponential::rv(rate, *rng);
//
//    
//    // rescale the times
//    for (size_t i = 0; i < n; i++ )
//    {
//        (*times)[i] = (*times)[i] *  origin / lastEvent;
//    }
//    
//    // finally sort the times
//    std::sort(times->begin(), times->end());
//	
//    return (*times)[0];
}


/**
 * Swap the parameters held by this distribution.
 *
 * 
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void DiversityDependentPureBirthProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == initialSpeciation) 
    {
        initialSpeciation = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == capacity) 
    {
        capacity = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else 
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }
    
}
