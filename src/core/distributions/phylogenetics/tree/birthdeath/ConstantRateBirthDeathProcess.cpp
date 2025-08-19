#include <cmath>
#include <iosfwd>
#include <vector>

#include "Clade.h"
#include "ConstantRateBirthDeathProcess.h"
#include "BirthDeathForwardSimulator.h"
#include "BirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

ConstantRateBirthDeathProcess::ConstantRateBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                                                     const TypedDagNode<double> *r, const TypedDagNode<double> *mp, const std::string& ss, const std::vector<Clade> &ic, const std::string &cdt,
                                                     const std::vector<Taxon> &tn) : BirthDeathProcess( ra, r, mp, ss, ic, cdt, tn, NULL ),
    speciation( s ),
    extinction( e )
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



void ConstantRateBirthDeathProcess::redrawValue( SimulationCondition sim_condition )
{

    if ( sim_condition == SimulationCondition::MCMC )
    {
        if ( starting_tree == NULL )
        {
            simulateTree(true);
        }
    }
    else if ( sim_condition == SimulationCondition::VALIDATION && condition == "nTaxa" )
    {
        simulateTree(false);
    }
    else if ( sim_condition == SimulationCondition::VALIDATION )
    {
        
        BirthDeathForwardSimulator simulator;
        
        size_t num_epochs = 1;
        std::vector< std::vector<double> > tmp = std::vector< std::vector<double> >( num_epochs, std::vector<double>(1,0) );

        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setBurstProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = extinction->getValue();
        simulator.setExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setMassExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 1.0;
        simulator.setSamplingProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = speciation->getValue();
        simulator.setSpeciationRate( tmp );
        
        simulator.setTimeline( std::vector<double>(1,0) );
        
        
        simulator.setRootCategoryProbabilities( std::vector<double>(1,1) );
        
        Tree *my_tree = simulator.simulateTreeConditionTime( getOriginAge(), BirthDeathForwardSimulator::SIM_CONDITION::ROOT);
        
        // store the new value
        delete value;
        value = my_tree;
        taxa = value->getTaxa();
    }
    else
    {
        throw RbException("Uknown condition for simulating tree in birth-death process.");
    }
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
