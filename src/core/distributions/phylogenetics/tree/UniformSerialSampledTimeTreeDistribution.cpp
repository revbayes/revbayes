#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "UniformSerialSampledTimeTreeDistribution.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;

UniformSerialSampledTimeTreeDistribution::UniformSerialSampledTimeTreeDistribution( const TypedDagNode<double> *a, const std::vector<Taxon> &n, const bool &ur) : TypedDistribution<Tree>( new Tree() ),
    start_age( a ),
    taxa( n ),
	has_root_age( ur )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( start_age );
    
    num_taxa = taxa.size();
    
    simulateTree();
}


UniformSerialSampledTimeTreeDistribution::~UniformSerialSampledTimeTreeDistribution()
{
    
}




/* Clone function */
UniformSerialSampledTimeTreeDistribution* UniformSerialSampledTimeTreeDistribution::clone( void ) const
{
    return new UniformSerialSampledTimeTreeDistribution( *this );
}

void UniformSerialSampledTimeTreeDistribution::buildSerialSampledRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &nodes, const std::vector<double> &ages)
{

	// Get the rng
	RandomNumberGenerator* rng = GLOBAL_RNG;

	// make a vector of active nodes, and a vector of serial nodes
	std::vector<TopologyNode*> active_nodes;
	std::vector<TopologyNode*> extinct_nodes;

	for(int i = 0; i < nodes.size(); ++i)
	{
		if ( nodes.at(i)->getAge() == 0.0 )
		{
			active_nodes.push_back( nodes.at(i) );
		}
		else
		{
			extinct_nodes.push_back( nodes.at(i) );
		}
	}

	// loop backward through ages
	size_t num_ages = ages.size();
	double current_time = 0.0;

	for(int i = num_ages - 1; i >= 0; i--)
	{
		// get the current time
		double current_time = ages[i];

		// check if any extinct nods become active
		size_t num_extinct = extinct_nodes.size();
		for(int j = num_extinct - 1; j >= 0; --j)
		{
			if ( extinct_nodes.at(j)->getAge() < current_time )
			{
				// add the extinct node to the active nodes list, remove it from the extinct nodes list
				active_nodes.push_back( extinct_nodes.at(j) );
				extinct_nodes.erase( extinct_nodes.begin() + long(j) );
			}
		}

		// randomly draw one child (arbitrarily called left) node from the list of active nodes
		size_t left = static_cast<size_t>( floor( rng->uniform01() * active_nodes.size() ) );
		TopologyNode* leftChild = active_nodes.at(left);

		// remove the randomly drawn node from the list
		active_nodes.erase( active_nodes.begin() + long(left) );

		// randomly draw one child (arbitrarily called left) node from the list of active nodes
		size_t right = static_cast<size_t>( floor( rng->uniform01() * active_nodes.size() ) );
		TopologyNode* rightChild = active_nodes.at(right);

		// remove the randomly drawn node from the list
		active_nodes.erase( active_nodes.begin() + long(right) );

		// make sure ages are consistent

		// add the parent
		TopologyNode* parent = new TopologyNode(i + num_taxa);
		parent->addChild(leftChild);
		parent->addChild(rightChild);
		leftChild->setParent(parent);
		rightChild->setParent(parent);
		parent->setAge( current_time );
		active_nodes.push_back(parent);
		nodes.push_back(parent);

	}

}


/* Compute probability */
double UniformSerialSampledTimeTreeDistribution::computeLnProbability( void )
{
    // Variable declarations and initialization
    double lnProb = 0.0;
    double age = 0.0;
    if ( has_root_age == true )
    {
    	age = start_age->getValue();
        // we need to check that the root age matches
        if ( age != value->getRoot().getAge() )
        {
            return RbConstants::Double::neginf;
        }
    }
    else
    {
    	age = start_age->getValue();
    }

    // check that the ages are in correct chronological order
    // i.e., no child is older than its parent
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        const TopologyNode &the_node = *(*it);
        // MRM: I'm disabling sampled ancestors for this distribution for now
        // because I don't know if they make sense!

        if ( the_node.isSampledAncestorTip() == true )
        {
        	return RbConstants::Double::neginf;
        }

        if ( the_node.isRoot() == false )
        {
            if ( (the_node.getAge() - (*it)->getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( (the_node.getAge() - (*it)->getParent().getAge()) > 1E-6 && the_node.isSampledAncestorTip() == true )
            {
                return RbConstants::Double::neginf;
            }
        }
        else if ( has_root_age == false )
        {
        	if ( the_node.getAge() > age )
        	{
        		return RbConstants::Double::neginf;
        	}
        }
    }

    // compute the probability of the coalescent events

    // get the coalescence times
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
    std::vector<double> serial_times;
    size_t num_taxa_at_present = value->getNumberOfTips();
    for (size_t i = 0; i < value->getNumberOfTips(); ++i)
    {
        const TopologyNode& n = value->getTipNode( i );
        double a = n.getAge();
        if ( a > 0.0 ) {
        	serial_times.push_back(a);
            --num_taxa_at_present;
        }
    }

    // compute the coalescent probabilities
    size_t current_num_taxa = num_taxa_at_present;
    double current_time = 0.0;
	for(int i = 0; i < ages.size(); ++i)
	{
		// get the current time
		double current_time = ages[i];

		// check if any extinct nodes become active
		size_t num_extinct = serial_times.size();

		for(int j = serial_times.size() - 1; j >= 0; --j)
		{
			if ( serial_times.at(j) < current_time )
			{
				// add the extinct node to the current number of taxa, remove it from the serial times list
				current_num_taxa++;
				serial_times.erase( serial_times.begin() + long(j) );
			}
		}

		// compute the probability of the coalescent event
		lnProb += std::log(2.0) - std::log( static_cast<double>(current_num_taxa) ) - std::log( static_cast<double>(current_num_taxa - 1) );

		// decrement the number of lineages after the coalescent event
		current_num_taxa--;
	}

	// compute the node-age probabilities
	size_t num_coalesced_nodes = value->getNumberOfTips() - 1;
	if ( has_root_age == false )
	{
		++num_coalesced_nodes;
	}
	for(size_t i = 1; i < num_coalesced_nodes; ++i)
	{
        const TopologyNode& n = value->getTipNode( i );
        double a = n.getAge();
        lnProb -= std::log(age - a);
	}

    return lnProb;
}


void UniformSerialSampledTimeTreeDistribution::redrawValue( void )
{
    simulateTree();
}

/**
 * Simulate new coalescent times.
 *
 * \param[in]    n      The number of coalescent events to simulate.
 *
 * \return    A vector of the simulated coalescent times.
 */
std::vector<double> UniformSerialSampledTimeTreeDistribution::simulateCoalescentAges() const
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // create the container for the ages
    std::vector<double> ages;

    // for each tip, simulate an age between max age and tip age
    double max_age = start_age->getValue();
    size_t num_ages;
    if ( has_root_age == true )
    {
    	// root age provided
    	ages.push_back(max_age);
    	num_ages = num_taxa - 2;
    }
    else
    {
    	// stem age provided
    	num_ages = num_taxa - 1;
    }

    for(size_t i = 0; i < num_ages; ++i)
    {
    	// get the age of the tip
    	double a = taxa[i + 1].getAge();

    	// simulate the age of a node
    	double new_age = a + rng->uniform01() * (max_age - a);

    	//add the age to the vector of ages
    	ages.push_back(new_age);
    }

    // sort the ages (from youngest to oldest)
    std::sort(ages.begin(), ages.end(), std::greater<double>());

    // return the ages
    return ages;
}

/** Simulate the tree conditioned on the time of origin */
void UniformSerialSampledTimeTreeDistribution::simulateTree( void )
{
    
    delete value;
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Create the time tree object (topology + times)
    Tree *psi = new Tree();

    // Root the topology by setting the appropriate flag
    psi->setRooted( true );
    
    // make a vector of tip nodes
    std::vector<TopologyNode* > nodes;

    // set tip names
    for (size_t i=0; i<num_taxa; ++i)
    {
        // get the node from the list
        TopologyNode* node = new TopologyNode(i);

        // set name and age
        const std::string& name = taxa[i].getName();
        node->setName(name);
        node->setSpeciesName(taxa[i].getSpeciesName());
        node->setAge(taxa[i].getAge());
        node->setTaxon(taxa[i]);
        // add to tips
        nodes.push_back(node);
    }

    // get times for simulation
    std::vector<double> ages = simulateCoalescentAges();

    // recursively build the tree
    buildSerialSampledRandomBinaryTree(psi, nodes, ages);

    // initialize the topology by setting the root
    TopologyNode* root = nodes[nodes.size() - 1];
    psi->setRoot(root, true);

    // finally store the new value
    value = psi;

}

void UniformSerialSampledTimeTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == start_age)
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}

/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void UniformSerialSampledTimeTreeDistribution::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == start_age )
    {
        dag_node->keepAffected();
    }
    
}

/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UniformSerialSampledTimeTreeDistribution::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == start_age )
    {
    	if ( has_root_age == true )
    	{
    		value->getNode( value->getRoot().getIndex() ).setAge( start_age->getValue() );
    	}
        dag_node->restoreAffected();
    }
    
}

/** Swap a parameter of the distribution */
void UniformSerialSampledTimeTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == start_age)
    {
    	start_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void UniformSerialSampledTimeTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == start_age )
    {
    	if ( has_root_age == true )
    	{
            value->getNode( value->getRoot().getIndex() ).setAge( start_age->getValue() );
    	}
        dag_node->touchAffected();
    }
    
}
