#include "Clade.h"
#include "ConstantRateBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

ConstantRateBirthDeathProcess::ConstantRateBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                                                     const TypedDagNode<double> *r, const std::string& ss, const std::vector<Clade> &ic, const std::string &cdt,
                                                     const std::vector<Taxon> &tn) : BirthDeathProcess( ra, r, ss, ic, cdt, tn ),
    speciation( s ),
    extinction( e )
{
    addParameter( speciation );
    addParameter( extinction );

    simulateTree();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ConstantRateBirthDeathProcess* ConstantRateBirthDeathProcess::clone( void ) const
{
    
    return new ConstantRateBirthDeathProcess( *this );
}



double ConstantRateBirthDeathProcess::lnSpeciationRate(double t) const
{
    double ln_lambda = log( speciation->getValue() );
    return ln_lambda;
}


double ConstantRateBirthDeathProcess::computeProbabilitySurvival(double start, double end) const 
{
    
    // compute the rate
    double mu = extinction->getValue();
    double lambda = speciation->getValue();
    double rate = mu - lambda;
    
    // do the integration of int_{t_low}^{t_high} ( mu(s) exp(rate(t,s)) ds )
    // where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx )
    
//    double den = 1.0 + ( exp(-rate*start) * mu / rate ) * ( exp(rate*end) - exp(rate*start) );
    double den = 1.0 + mu / rate * ( exp(rate*(end-start)) - 1 );
    
    return (1.0 / den);
}



double ConstantRateBirthDeathProcess::rateIntegral(double t_low, double t_high) const
{
    
    double rate = (speciation->getValue() - extinction->getValue()) * (t_low - t_high);
        
    return rate;
}



double ConstantRateBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the parameters
    double age = origin - present;
    double b = speciation->getValue();
    double d = extinction->getValue();
    double r = rho->getValue();
 
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
void ConstantRateBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
        BirthDeathProcess::swapParameterInternal(oldP, newP);
    }
    
}
