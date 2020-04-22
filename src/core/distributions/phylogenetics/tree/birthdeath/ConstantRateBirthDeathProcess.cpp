#include <cmath>
#include <iosfwd>
#include <vector>

#include "ConstantRateBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "DistributionPoisson.h"
#include "DistributionExponential.h"
#include "BirthDeathProcess.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }

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
    double f = rho->getValue();
    double rate = mu - lambda;
    
    // do the integration of int_{t_low}^{t_high} ( mu(s) exp(rate(t,s)) ds )
    // where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx )
    
//    double den = 1.0 + ( exp(-rate*start) * mu / rate ) * ( exp(rate*end) - exp(rate*start) );
    double den = 1.0 + mu / rate * ( exp(rate*(end-start)) - 1 );
    
    return (1.0 / den);
}

double ConstantRateBirthDeathProcess::computeProbabilityNoSurvival(double start, double end) const
{
    double mu = extinction->getValue();
    double lambda = speciation->getValue();
    double f = rho->getValue();
    //    double rate = mu - lambda;
    
    double p = 1.0 - f*(lambda-mu) / (f*lambda + (lambda*(1.0-f)-mu) * std::exp(-(lambda-mu)*(start-end)));
    
    return p;
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


bool ConstantRateBirthDeathProcess::simulateOneHiddenCladeRecursively(TopologyNode* curr_node, std::vector<TopologyNode*>& stored_nodes)
{
    
    bool success = true;
    
    // model parameters
    double lambda = speciation->getValue();
    double mu = extinction->getValue();
    double total_rate = lambda + mu;
    double sampling_prob = rho->getValue();
    
    // sample the next event
    double dt = RbStatistics::Exponential::rv(total_rate, *GLOBAL_RNG);
    double next_age = curr_node->getAge() - dt;
    
    
    // does the lineage survive until present?
    bool is_sampled = (GLOBAL_RNG->uniform01() < sampling_prob);
    if (next_age < 0.0)
    {
        if (is_sampled) {
//            std::cout << "        SAMPLED " << 0.0 << "\n";
            return false;
        }
        else {
//            std::cout << "        UNSAMPLED " << 0.0 << "\n";
            // next node
            std::stringstream ss_name;
            ss_name << "U" << (long)(GLOBAL_RNG->uniform01()*1e9);
            TopologyNode* next_node = new TopologyNode( ss_name.str() );
            next_node->setSampledDescendant(false);
            next_node->setAge(0.0);
            next_node->setParent(curr_node);
            curr_node->addChild(next_node);
            
            stored_nodes.push_back(next_node);
            return true;
        }
    }
    
    // next node
    std::stringstream ss_name;
    ss_name << "X" << (long)(GLOBAL_RNG->uniform01()*1e9);
    TopologyNode* next_node = new TopologyNode( ss_name.str() );
    next_node->setSampledDescendant(false);
    next_node->setAge(next_age);
    next_node->setParent(curr_node);
    curr_node->addChild(next_node);
    
    stored_nodes.push_back(next_node);
    
    // sample next event type
    bool is_speciation = (GLOBAL_RNG->uniform01() < (lambda/total_rate));
    if (is_speciation) {
//        std::cout << "        SPECIATION " << next_node->getAge() << "\n";
        success &= simulateOneHiddenCladeRecursively(next_node, stored_nodes);
        success &= simulateOneHiddenCladeRecursively(next_node, stored_nodes);
    } else {
        
//        std::cout << "        EXTINCTION " << next_node->getAge() << "\n";
        // do nothing for extinction
    }

    return success;
}

TopologyNode* ConstantRateBirthDeathProcess::simulateOneHiddenClade(double time) {

    
    TopologyNode* new_clade = NULL;
    int num_attempts = 10000;
    int count = 1;
    bool success = false;
    
//    std::cout << "  HIDDEN_CLADE " << time << "\n";
    while (!success && count <= num_attempts) {

//        std::cout << "    BRANCH-ATTEMPT " << count << "\n";
        std::vector<TopologyNode*> stored_nodes;
        std::stringstream ss_name;
        ss_name << "node_" << (long)(GLOBAL_RNG->uniform01()*1e16);
        new_clade = new TopologyNode( ss_name.str() );
        new_clade->setAge(time);
        new_clade->setSampledDescendant(false);
        stored_nodes.push_back(new_clade);
        
        success = simulateOneHiddenCladeRecursively(new_clade, stored_nodes);
        
        new_clade->setTree( this->value );
        
        // delete all nodes if simulation fails
        if (!success) {
            delete new_clade;
            stored_nodes.clear();
        }
        
        count += 1;
    }
    
    if (success) {
        return new_clade;
    } else {
        return NULL;
    }
}

/**
 * Simulate hidden clades upon the reconstructed tree under a birth-death process.
 **/

void ConstantRateBirthDeathProcess::simulateHiddenClades(void) {
    
    bool success = false;
    int max_attempts = 100;
    int num_attempts = 0;
    
    
    Tree restore_tree( *this->value );
    while (!success && num_attempts <= max_attempts) {
    
        // update tree
        this->value->pruneTaxaWithoutSampledDescendants();
//        this->value = &restore_tree;
        
        // simulate hidden speciation events along reconstructed tree
        double birth_rate = speciation->getValue();
        double tree_length = value->getTreeLength();
        double max_time = value->getRoot().getAge();
        size_t num_events = RbStatistics::Poisson::rv(birth_rate * tree_length, *GLOBAL_RNG);
        size_t num_sampled_events = 0;
//        std::cout << "TREE-ATTEMPT " << num_attempts << "\n";
        
        // simulate times for each hidden speciation event
        std::vector<TopologyNode*> nodes = value->getNodes();
        std::map<TopologyNode*, std::vector<double> > thinned_events;
        
        // attempt to place each event
        for (size_t i = 0; i < num_events; i++)
        {
            // sample a branch from the tree w.p. branch_length / tree_length
            double t = tree_length * GLOBAL_RNG->uniform01();
            size_t node_index = 0;
            while (t > 0.0)
            {
                // decrement by branch lengths until t is negative
                TopologyNode* nd = nodes[node_index];
                double branch_length = nd->getBranchLength();
                t -= branch_length;
                
                // have we found the branch?
                if (t <= 0.0 ) {
                    
                    // once the branch is sampled, select a Poisson random time
                    double event_time = branch_length * GLOBAL_RNG->uniform01();
                    
                    // only add the event if it leaves no extant samples
                    double p_no_extant_samples = computeProbabilityNoSurvival(-event_time, 0.0); //max_time);
                    bool has_no_extant_samples = (GLOBAL_RNG->uniform01() < p_no_extant_samples);
                    if (has_no_extant_samples) {
                        thinned_events[nd].push_back( event_time );
                        num_sampled_events += 1;
                    }
                    
                }
                else {
                    // try next branch
                    node_index++;
                }
            }
        }
//        std::cout << "  num_total_events = " << num_events << "\n";
//        std::cout << "  num_sampled_events = " << num_sampled_events << "\n";
        
        // MJL this may be easier to do when sampling events & times
        
        // now simulate subclades
        bool success_subclade = true;
        std::map< TopologyNode*, std::vector<double> >::iterator it;
        for (it = thinned_events.begin(); it != thinned_events.end(); it++)
        {
            TopologyNode* end_node = it->first;
            TopologyNode* start_node = &end_node->getParent();
            std::vector<double> times = it->second;
            std::sort( times.begin(), times.end() );
            
            for (int i = (int)times.size() - 1; i >= 0; i--) {
                
                // generate subclade
                TopologyNode* new_subclade = simulateOneHiddenClade( end_node->getAge() + times[i] );
                
                // did it succeed?
                success_subclade &= (new_subclade != NULL);
                
                // new attempt upon failure
                if (!success_subclade) break;
                
                // otherwise graft the subclade
                new_subclade->setParent( start_node );
                new_subclade->addChild( end_node );
                start_node->removeChild( end_node );
                start_node->addChild( new_subclade );
                end_node->setParent( new_subclade );
                start_node = new_subclade;
                
            }
            // new tree-attempt upon failure
            if (!success_subclade) {
                success = false;
                break;
            }
            else {
                success = true;
            }
        }
        
        num_attempts += 1;
    }
    
    if (success) {
//        std::vector<TopologyNode*> nodes = value->getNodes();
//        for (size_t i = 0; i < taxa.size(); i++) {
//            TopologyNode* n = nodes[i];
//            if (n->isTip()) {
//                std::vector<Taxon>::iterator it = std::find( taxa.begin(), taxa.end(), n->getTaxon() );
//                if (it == taxa.end()) {
//                    taxa.
//                }
//            }
//        }
        value->setRoot( &value->getRoot(), true );
        value->getNewickRepresentation();
        value->clearTaxonBitSetMap();
        value->getTaxonBitSetMap();
//        std::cout << value << "\n";
//        std::cout << *value << "\n";
    }
    else {
        throw RbException("Failed to generate hidden clades after 10000 attempts.");
    }
    
    return;
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
