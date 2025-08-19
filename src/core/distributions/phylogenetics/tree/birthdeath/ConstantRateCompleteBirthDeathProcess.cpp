#include <cstddef>
#include <cmath>
#include <iosfwd>
#include <vector>

#include "ConstantRateCompleteBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
#include "AbstractBirthDeathProcess.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

ConstantRateCompleteBirthDeathProcess::ConstantRateCompleteBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<double> *e, const std::string &cdt, const std::vector<Taxon> &tn, bool uo, bool mr) : AbstractBirthDeathProcess( ra, cdt, tn, uo, NULL ),
    speciation( s ),
    extinction( e ),
    use_mrca( mr )
{
    addParameter( speciation );
    addParameter( extinction );

    simulateTree(true);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ConstantRateCompleteBirthDeathProcess* ConstantRateCompleteBirthDeathProcess::clone( void ) const
{

    return new ConstantRateCompleteBirthDeathProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return   The log-transformed probability density.
 */
double ConstantRateCompleteBirthDeathProcess::computeLnProbabilityTimes( void ) const
{

    double ln_prob_times = 0;
    size_t num_extant_left  = num_extant(value->getRoot().getChild(0), &ln_prob_times);
    size_t num_extant_right = num_extant(value->getRoot().getChild(1), &ln_prob_times);

    if ( use_mrca == true && (num_extant_left == 0 || num_extant_right == 0) )
    {
        return RbConstants::Double::neginf;
    }

    if ( use_mrca == false )
    {
        ln_prob_times += getOriginAge() - getRootAge();
    }

    ln_prob_times *= - ( speciation->getValue() + extinction->getValue() );

    ln_prob_times += ( value->getNumberOfTips() - 1 ) * log( speciation->getValue() );


    size_t num_extinct = value->getNumberOfTips() - num_extant_left - num_extant_right;

    if ( extinction->getValue() > 0.0 )
    {
        ln_prob_times += num_extinct * log( extinction->getValue() );
    }

    return ln_prob_times;
}


size_t ConstantRateCompleteBirthDeathProcess::num_extant( const TopologyNode & node, double* tl ) const
{
    size_t num = 0;

    if ( node.isRoot() == false )
    {
        *tl += node.getParent().getAge() - node.getAge();
    }

    if (node.isTip())
    {
        num = (node.getAge() == 0);
    }
    else
    {
        for (size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            num += num_extant(node.getChild(i), tl);
        }
    }

    return num;
}

/**
 *
 */
double ConstantRateCompleteBirthDeathProcess::lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const
{
    for ( size_t i = 0; i < n; i++)
    {
        if ( value->getNode(i).getAge() > 0 )
        {
            n--;
        }
    }

    double b = speciation->getValue();
    double d = extinction->getValue();
    double t = end - start;

    double e = exp(-(b - d)*t);

    double u = b*(1.0 - d*e)/(b - d*e);

    if ( n == 0 )
    {
        return log(1.0 - pSurvival(start,end));
    }

    double lnN = log(pSurvival(start,end)) + log(1.0 - u) + (n - 1)*log(u);

    return lnN;
}


double ConstantRateCompleteBirthDeathProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to the probability of a (ranked) non-oriented tree with labeled extant tips
    // by multiplying by 2^{n-1} / n!
    // where n is the number of tips

    int n = (int)value->getNumberOfTips();
    int m = (int)value->getNumberOfExtinctTips();

    // condition on survival
    if (n == m)
    {
        return RbConstants::Double::neginf;
    }

    return (n - 1) * RbConstants::LN2 - RbMath::lnFactorial(n - m);
}


/**
 * Compute the probabililty of survival (no extinction) of the process.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return The probability of survival of the process.
 */
double ConstantRateCompleteBirthDeathProcess::pSurvival(double start, double end) const
{
    double b = speciation->getValue();
    double d = extinction->getValue();
    double t = end - start;

    double pS = (b - d) / (b - d * exp(-(b - d)*t));

    return pS;
}


double ConstantRateCompleteBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the parameters
    double age = origin - present;
    double b = speciation->getValue();
    double d = extinction->getValue();
    double r = 1.0;
 
    // get a random draw
    double u = rng->uniform01();

    
    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    else
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    
    return present + t;
}


/** Swap a parameter of the distribution */
void ConstantRateCompleteBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == speciation)
    {
        speciation = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == extinction)
    {
        extinction = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }

}
