#include "BirthDeathBurstProcess.h"

#include <cmath>
#include <string>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
#include "AbstractRootedTreeDistribution.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }


using namespace RevBayesCore;


/**
 * Constructor.
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    ra             Age or the root (=time of the process).
 * \param[in]    s              Speciation rate.
 * \param[in]    e              Extinction rate.
 * \param[in]    p              Extinction sampling rate.
 * \param[in]    r              Sampling probability at present time.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Clades conditioned to be present.
 */
BirthDeathBurstProcess::BirthDeathBurstProcess( const TypedDagNode<double> *ra,
                                                const TypedDagNode<double> *s,
                                                const TypedDagNode<double> *e,
                                                const TypedDagNode<double> *b,
                                                const TypedDagNode<double> *bt,
                                                const TypedDagNode<double> *r,
                                                const std::string& cdt,
                                                const std::vector<Taxon> &tn) : AbstractBirthDeathProcess( ra, cdt, tn, false, NULL ),
    lambda( s ),
    mu( e ),
    beta( b ),
    time_burst( bt ),
    rho( r ),
    lineage_bursted_at_event(tn.size()*2-1,false)
{
    addParameter( lambda );
    addParameter( mu );
    addParameter( beta );
    addParameter( time_burst );
    addParameter( rho );
    
    simulateTree(true);
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
BirthDeathBurstProcess* BirthDeathBurstProcess::clone( void ) const
{
    
    return new BirthDeathBurstProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double BirthDeathBurstProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    // prepare the probability computation
    prepareProbComputation();
    
    // variable declarations and initialization
    double lnProbTimes = computeLnProbabilityTimes();
    
    return lnProbTimes;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return    The log-probability density.
 */
double BirthDeathBurstProcess::computeLnProbabilityTimes( void ) const
{
    
    double lnProbTimes = 0.0;
    double process_time = getOriginAge();
    size_t num_initial_lineages = 2;
    
    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double burst_prob = beta->getValue();
    double burst_time = time_burst->getValue();
    double sampling_prob = rho->getValue();
    
    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();
    size_t num_lineages_burst_at_event = 0;
    size_t num_lineages_alive_at_burst = 0;
    
    // classify nodes
    int num_extant_taxa = 0;
    
    std::vector<double> internal_node_ages = std::vector<double>();
    for (size_t i = 0; i < num_nodes; ++i)
    {
        const TopologyNode& n = value->getNode( i );
        double node_age = n.getAge();
        
        // store the age if this node is an internal node
        if ( n.isInternal() == true && n.isRoot() == false )
        {
            internal_node_ages.push_back( node_age );
        }
        else if ( n.isTip() == true )
        {
            ++num_extant_taxa;
        }
        
        // count the number of lineages alive at the burst event
        // do not count the lineages that actually bursted
        if ( n.isRoot() == false )
        {
            double parent_age = n.getParent().getAge();
            if ( fabs(node_age - burst_time) > 1E-4 &&
                 node_age < burst_time &&
                 fabs(parent_age - burst_time) > 1E-4 &&
                 parent_age > burst_time )
            {
                ++num_lineages_alive_at_burst;
            }
            
        }
        
        // count the number of lineages that bursted
        if ( lineage_bursted_at_event[i] == true )
        {
            
            if ( fabs(node_age - burst_time) < 1E-4 )
            {
                // node is sampled ancestor
                ++num_lineages_burst_at_event;
            }
            else
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
        
    }
    
    // add the log probability for sampling the extant taxa
    lnProbTimes += num_extant_taxa * log( sampling_prob );
    
    // add the log probability of the initial sequences
    lnProbTimes += lnQ(process_time) * num_initial_lineages;
    
    // add the log probability for the burst event
    if ( num_lineages_burst_at_event > 0 )
    {
        lnProbTimes += log(burst_prob) * num_lineages_burst_at_event;
        lnProbTimes += log(1.0-burst_prob+2*burst_prob*pZero(burst_time)) * num_lineages_alive_at_burst;
    }
    
    // add the log probability for the internal node ages
    lnProbTimes += (internal_node_ages.size()-num_lineages_burst_at_event) * log( birth_rate );
    for (size_t i=0; i<internal_node_ages.size(); i++)
    {
        lnProbTimes += lnQ(internal_node_ages[i]);
    }
    
    // condition on survival
    if ( condition == "survival")
    {
        lnProbTimes -= num_initial_lineages * log(1.0 - pZero(process_time));
    }
    // condition on nTaxa
    else if ( condition == "nTaxa" )
    {
        lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
    }
    
    
    return lnProbTimes;
    
}


void BirthDeathBurstProcess::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, long &rv) const
{
    
    if ( name == "numBurstSpeciations" )
    {
        rv = 0;
        size_t num_nodes = value->getNumberOfNodes();
        
        for (size_t i=0; i<num_nodes; ++i)
        {
            
            if ( lineage_bursted_at_event[i] == true )
            {
                ++rv;
            }
            
        }
        
    }
    else
    {
        throw RbException() << "The character dependent birth-death process does not have a member method called '" <<  name << "'.";
    }
    
}


/**
 * Get the affected nodes by a change of this node.
 * If the burst age has changed than we need to call get affected again.
 */
void BirthDeathBurstProcess::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == time_burst )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
    // delegate to base class
    AbstractRootedTreeDistribution::getAffected(affected, affecter);
    
}


double BirthDeathBurstProcess::getBurstTime( void ) const
{
    
    return time_burst->getValue();
}


bool BirthDeathBurstProcess::isBurstSpeciation( size_t index ) const
{
    
    return lineage_bursted_at_event[index];
}


/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void BirthDeathBurstProcess::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == time_burst && dag_node != NULL)
    {
        dag_node->keepAffected();
    }
    
    // delegate to base class
    AbstractRootedTreeDistribution::keepSpecialization( affecter );
    
}


double BirthDeathBurstProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n-1} / n!
    // where n is the number of extant tips
    
    int num_taxa = (int)value->getNumberOfTips();
    
    return (num_taxa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa);
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Speciation rate at time t.
 */
double BirthDeathBurstProcess::pSurvival(double start, double end) const
{
    return 1.0 - pZero(end);
}


/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void BirthDeathBurstProcess::restoreSpecialization(const DagNode *affecter)
{
    if ( affecter == time_burst )
    {
        
        size_t num_nodes = value->getNumberOfNodes();
        double new_burst_time = time_burst->getValue();
        
        for (size_t i=0; i<num_nodes; ++i)
        {
            
            if ( lineage_bursted_at_event[i] == true )
            {
                value->getNode(i).setAge( new_burst_time );
            }
        }
        
        if ( dag_node != NULL )
        {
            dag_node->restoreAffected();
        }
    }
    
    // delegate to base class
    AbstractRootedTreeDistribution::restoreSpecialization( affecter );
    
}


void BirthDeathBurstProcess::setBurstSpeciation( size_t index, bool tf )
{
    
    lineage_bursted_at_event[index] = tf;
}



void BirthDeathBurstProcess::setValue(Tree *v, bool force)
{
    AbstractBirthDeathProcess::setValue(v,force);
    
    double t = getBurstTime();
    
    lineage_bursted_at_event.resize( v->getNumberOfNodes() );
    
    for (size_t i=0; i<value->getNumberOfNodes(); ++i)
    {
        const TopologyNode& n = value->getNode(i);
        if ( n.isTip() == false )
        {
            setBurstSpeciation(i, fabs( n.getAge() - t ) < 1E-4);
        }
    }
}


/**
 * Simulate new speciation times.
 */
double BirthDeathBurstProcess::simulateDivergenceTime(double origin, double present) const
{
    
    // incorrect placeholder
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the parameters
    double age = origin - present;
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho->getValue();
    
    // get a random draw
    double u = rng->uniform01();
    
    
    // compute the time for this draw
    double t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    
    return present + t;
}


double BirthDeathBurstProcess::pZero(double t) const
{
    double birth = lambda->getValue();
    double death = mu->getValue();
    double burst = beta->getValue();
    double t_b = time_burst->getValue();
    double sampling = rho->getValue();
    double E = 0;
    
    double A = birth - death;
    
    if ( t < (t_b+1E-4) )
    {
        double B = ((1.0 - 2.0*(1.0-sampling)) * birth + death ) /  A;
        E = birth + death - A * (1.0 + B - exp(-A*t) * (1.0-B)) / (1.0+B+exp(-A*t)*(1.0-B));
        E /= (2*birth);
    }
    else
    {
        double B_tmp = ((1.0 - 2.0*(1.0-sampling)) * birth + death ) /  A;
        double E_tmp = birth + death - A * (1.0 + B_tmp - exp(-A*t_b) * (1.0-B_tmp)) / (1.0+B_tmp+exp(-A*t_b)*(1.0-B_tmp));
        E_tmp /= (2*birth);
        double B = ((1.0 - 2.0*((1.0-burst)*E_tmp+burst*E_tmp*E_tmp)) * birth + death ) /  A;
        E = birth + death - A * (1.0 + B - exp(-A*(t-t_b)) * (1.0-B)) / (1.0+B+exp(-A*(t-t_b))*(1.0-B));
        E /= (2*birth);
    }
    
    return E;
}


double BirthDeathBurstProcess::lnQ(double t) const
{
    // D(t) = \frac{ 4e^{-A(t-t_y)} }{ \big(1+B_{t_y}+e^{-A(t-t_y)}(1-B_{t_y}) \big)^2 }
    
    double birth = lambda->getValue();
    double death = mu->getValue();
    double burst = beta->getValue();
    double t_b = time_burst->getValue();
    double sampling = rho->getValue();
    
    double A = birth - death;
    double B = 0.0;
    double D = 0.0;
    if ( t < (t_b+1E-4) )
    {
        B = ((1.0 - 2.0*(1.0-sampling)) * birth + death ) /  A;
        
        D = 4.0*exp(-A*t);
        double tmp = 1.0 + B + exp(-A*t)*(1.0-B);
        D /= (tmp*tmp);
    }
    else
    {
        double B_previous = ((1.0 - 2.0*(1.0-sampling)) * birth + death ) /  A;
        double E_previous = birth + death - A * (1.0 + B_previous - exp(-A*t_b) * (1.0-B_previous)) / (1.0+B_previous+exp(-A*t_b)*(1.0-B_previous));
        E_previous /= (2*birth);

        B = ((1.0 - 2.0*((1.0-burst)*E_previous+burst*E_previous*E_previous)) * birth + death ) /  A;

        
        double lnD_previous = 2*RbConstants::LN2 + (-A*t_b);
        lnD_previous -= 2 * log(1 + B_previous + exp(-A * t_b) * (1 - B_previous));
        
        double this_lnD = lnD_previous;
        this_lnD += 2*RbConstants::LN2 + (-A * (t - t_b));
        this_lnD -= 2 * log(1 + B + exp(-A * (t - t_b)) * (1 - B));

        return this_lnD;
    }
    
    return log( D );
}


#include "BirthDeathForwardSimulator.h"
void BirthDeathBurstProcess::redrawValue( SimulationCondition condition )
{

    if ( condition == SimulationCondition::MCMC )
    {
        if ( starting_tree == NULL )
        {
            simulateTree(true);
        }
    }
    else if ( condition == SimulationCondition::VALIDATION )
    {
        
        BirthDeathForwardSimulator simulator;
        
        size_t num_epochs = 2;
        std::vector< std::vector<double> > tmp = std::vector< std::vector<double> >( num_epochs, std::vector<double>(1,0) );

        tmp[0][0] = 0.0;
        tmp[1][0] = beta->getValue();
        simulator.setBurstProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = mu->getValue();
        simulator.setExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setMassExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        tmp[0][0] = 1.0;
        simulator.setSamplingProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setSamplingExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = lambda->getValue();
        simulator.setSpeciationRate( tmp );
        
        std::vector<double> timeline = std::vector<double>(2,0);
        timeline[1] = getBurstTime();
        simulator.setTimeline( timeline );
        
        
        
        simulator.setRootCategoryProbabilities( std::vector<double>(1,1) );
        
        Tree *my_tree = simulator.simulateTreeConditionTime( getOriginAge(), BirthDeathForwardSimulator::SIM_CONDITION::ROOT);
        
        // store the new value
//        delete value;
        setValue(my_tree);
        taxa = value->getTaxa();
    }
    else
    {
        throw RbException("Uknown condition for simulating tree in birth-death-burst process.");
    }
}



/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void BirthDeathBurstProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == time_burst)
    {
        time_burst = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    // delegate the super-class
    AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    
}



/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void BirthDeathBurstProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == time_burst )
    {
        
        size_t num_nodes = value->getNumberOfNodes();
        double new_burst_time = time_burst->getValue();
        
        for (size_t i=0; i<num_nodes; ++i)
        {
            
            if ( lineage_bursted_at_event[i] == true )
            {
                value->getNode(i).setAge( new_burst_time );
            }
        }
        
        if ( dag_node != NULL )
        {
            dag_node->touchAffected();
        }
    }
    
    // delegate to base class
    AbstractRootedTreeDistribution::touchSpecialization( affecter, touchAll );
    
}


