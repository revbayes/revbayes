#include <cstddef>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "AbstractMultispeciesCoalescent.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
#include "RbException.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

AbstractMultispeciesCoalescent::AbstractMultispeciesCoalescent(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : TypedDistribution<Tree>( NULL ),
    taxa(t),
    species_tree( sp ),
    num_taxa( taxa.size() ),
    log_tree_topology_prob (0.0)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( species_tree );

    std::set<std::string> species_names;
    for (std::vector<Taxon>::const_iterator it=taxa.begin(); it!=taxa.end(); ++it)
    {
        species_names.insert( it->getSpeciesName() );
    }

    double ln_fact = RbMath::lnFactorial((int)(num_taxa));

    log_tree_topology_prob = (num_taxa - 1) * RbConstants::LN2 - 2.0 * ln_fact - std::log( num_taxa ) ;

    redrawValue();

}


AbstractMultispeciesCoalescent::~AbstractMultispeciesCoalescent()
{

}



void AbstractMultispeciesCoalescent::attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times) {

    if (index < num_taxa-1)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;

        // randomly draw one node from the list of tips
        size_t tip_index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );

        // get the node from the list
        TopologyNode* parent = tips.at(tip_index);
        psi->getNode( parent->getIndex() ).setAge( times[num_taxa - index - 2] );

        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+tip_index);

        // add a left child
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            tips.push_back(leftChild);
        }

        // add a right child
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            tips.push_back(rightChild);
        }

        // recursive call to this function
        attachTimes(psi, tips, index+1, times);
    }
}


void AbstractMultispeciesCoalescent::buildRandomBinaryTree(std::vector<TopologyNode*> &tips)
{

    if (tips.size() < num_taxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;

        // randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );

        // get the node from the list
        TopologyNode* parent = tips.at(index);

        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+index);

        // add a left child
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);

        // add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);

        // recursive call to this function
        buildRandomBinaryTree(tips);
    }
}


double AbstractMultispeciesCoalescent::computeLnProbability( void )
{
    resetTipAllocations();

    // variable declarations and initialization
    double ln_prob_coal = 0;

    const Tree &sp = species_tree->getValue();

    ln_prob_coal = recursivelyComputeLnProbability( sp.getRoot() );

    return ln_prob_coal; // + logTreeTopologyProb;

}


double AbstractMultispeciesCoalescent::drawNe( size_t index )
{

    return 1.0;
}


double AbstractMultispeciesCoalescent::recursivelyComputeLnProbability( const RevBayesCore::TopologyNode &species_node )
{

    double ln_prob_coal = 0;

    if ( species_node.isTip() == false )
    {
        individuals_per_branch[ species_node.getIndex() ].clear();

        for (size_t i=0; i<species_node.getNumberOfChildren(); ++i)
        {
            ln_prob_coal += recursivelyComputeLnProbability( species_node.getChild(i) );
        }
    }

//    const TopologyNode *species_parent_node = NULL;
    double species_age = species_node.getAge();
    double parent_species_age = RbConstants::Double::inf;

    //        double branchLength = RbConstants::Double::inf;
    if ( species_node.isRoot() == false )
    {
        const TopologyNode &species_parent_node = species_node.getParent();
        parent_species_age = species_parent_node.getAge();
    }

    // create a local copy of the individuals per branch
    const std::set<const TopologyNode*> &initial_individuals = individuals_per_branch[species_node.getIndex()];
    std::set<const TopologyNode*> remaining_individuals = initial_individuals;

    // get all coalescent events among the individuals
    std::vector<double> coal_times;


    std::map<double, const TopologyNode *> coal_times_2_nodes;
    for ( std::set<const TopologyNode*>::iterator it = remaining_individuals.begin(); it != remaining_individuals.end(); ++it)
    {
        const TopologyNode *ind = (*it);
        if ( ind->isRoot() == false )
        {
            const TopologyNode &parent = ind->getParent();
            double parent_age = parent.getAge();
            coal_times_2_nodes[ parent_age ] = &parent;
        }
    }

    double current_time = species_age;
    while ( current_time < parent_species_age && coal_times_2_nodes.size() > 0 )
    {

        const TopologyNode *parent = coal_times_2_nodes.begin()->second;
        double parent_age = parent->getAge();
        current_time = parent_age;

        if ( parent_age < parent_species_age )
        { //Coalescence in the species tree branch

            // get the left and right child of the parent
            const TopologyNode *left = &parent->getChild( 0 );
            const TopologyNode *right = &parent->getChild( 1 );
            if ( remaining_individuals.find( left ) == remaining_individuals.end() || remaining_individuals.find( right ) == remaining_individuals.end() )
            {
                // one of the children does not belong to this species tree branch
                return RbConstants::Double::neginf;
            }

            //We remove the coalescent event and the coalesced lineages
            coal_times_2_nodes.erase( coal_times_2_nodes.begin() );
            remaining_individuals.erase( remaining_individuals.find( left ) );
            remaining_individuals.erase( remaining_individuals.find( right ) );

            //We insert the parent in the vector of lineages in this branch
            remaining_individuals.insert( parent );
            if ( parent->isRoot() == false )
            {
                const TopologyNode *grand_parent = &parent->getParent();
                coal_times_2_nodes[ grand_parent->getAge() ] = grand_parent;
            }

            coal_times.push_back( parent_age );


        } //End if coalescence in the species tree branch
        else
        { //No more coalescences in this species tree branch

            // jump out of the while loop
            //                currentTime = speciesAge;
            break;
        }


    } // end of while loop

    if ( initial_individuals.size() > 1 )
    {
        ln_prob_coal += computeLnCoalescentProbability(initial_individuals.size(), coal_times, species_age, parent_species_age, species_node.getIndex(), species_node.isRoot() == false);
    }


    // merge the two sets of individuals that go into the next species
    if ( species_node.isRoot() == false )
    {
        std::set<const TopologyNode *> &incoming_lineages = individuals_per_branch[ species_node.getParent().getIndex() ];
        incoming_lineages.insert( remaining_individuals.begin(), remaining_individuals.end());
    }

    return ln_prob_coal;
}


void AbstractMultispeciesCoalescent::redrawValue( void )
{

    simulateTree();

}


void AbstractMultispeciesCoalescent::resetTipAllocations( void )
{

    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();

    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (auto& node: species_tree_nodes)
    {
        if ( node->isTip() )
        {
            const std::string &name = node->getName();
            species_names_2_species_nodes[name] = node;
        }
    }

    // second, let's create a map from individual names to the species names
    std::map<std::string, std::string> individual_names_2_species_names;
    for (auto& taxon: taxa)
        individual_names_2_species_names[taxon.getName()] = taxon.getSpeciesName();

    // create a map for the individuals to branches
    individuals_per_branch = std::vector< std::set< const TopologyNode* > >(sp.getNumberOfNodes() );
    for (size_t i=0; i<num_taxa; ++i)
    {
//        const std::string &tip_name = it->getName();
        const TopologyNode &n = value->getNode( i );
        const std::string &individual_name = n.getName();
        const std::string &species_name = individual_names_2_species_names.at(individual_name);

        TopologyNode *species_node;
        try
        {
            species_node = species_names_2_species_nodes.at(species_name);
        }
        catch (...)
        {
            throw RbException()<<"No tree node for species '"<<species_name<<"'";
        }
        individuals_per_branch[ species_node->getIndex() ].insert( &n );
    }


}


/**
 * Set the current value.
 */
void AbstractMultispeciesCoalescent::setValue(Tree *v, bool f )
{

    // delegate to super class
    TypedDistribution<Tree>::setValue(v, f);

    resetTipAllocations();

}




void AbstractMultispeciesCoalescent::simulateTree( void )
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            species_names_2_nodes[name] = *it;
        }
    }


    std::map< const TopologyNode *, std::vector< TopologyNode* > > individuals_per_branch;

    for (std::vector< Taxon >::iterator it = taxa.begin(); it != taxa.end(); ++it)
    {
        TopologyNode *n = new TopologyNode( *it );
        const std::string &species_name = n->getSpeciesName();

        if ( species_name == "" )
        {
            throw RbException("Cannot match a taxon without species to a tip in the species tree. The taxon map is probably wrong.");
        }

        TopologyNode *species_node = species_names_2_nodes.at(species_name);

        if ( species_node == NULL )
        {
            throw RbException() << "Could not match a taxon with name" << species_name << " to any of the tips in the species tree.";
        }

        n->setAge( 0.0 );
        std::vector< TopologyNode * > &nodes_at_node = individuals_per_branch[ species_node ];
        nodes_at_node.push_back( n );
    }


    /*   for (std::map<std::string,std::string>::iterator it = gene2species.begin(); it != gene2species.end(); ++it)
     {
     const std::string &tipName = it->first;
     TopologyNode *n = new TopologyNode( tipName );
     const std::string &speciesName = it->second;
     TopologyNode *speciesNode = speciesNames2Nodes[speciesName];
     std::vector< TopologyNode * > &nodesAtNode = individualsPerBranch[ speciesNode ];
     nodesAtNode.push_back( n );
     }
     */

    std::map<TopologyNode *, double> nodes_2_ages;
    TopologyNode *root = NULL;
    // we loop over the nodes of the species tree in phylogenetic traversal
    for (std::vector<TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        TopologyNode *sp_node = *it;
        const TopologyNode *sp_parent_node = NULL;
        double branch_length = RbConstants::Double::inf;
        if ( sp_node->isRoot() == false )
        {
            sp_parent_node = &sp_node->getParent();
            branch_length = sp_parent_node->getAge() - sp_node->getAge();
        }

        std::vector<TopologyNode*> initial_individuals_at_branch = individuals_per_branch[sp_node];
        double branch_ne = drawNe( sp_node->getIndex() );

        double theta = 2.0 / branch_ne;

        double prev_coalescent_time = 0.0;

        size_t j = initial_individuals_at_branch.size();
        double n_pairs = j * (j-1) / 2.0;
        double lambda = n_pairs * theta;
        double u = RbStatistics::Exponential::rv( lambda, *rng);
        double next_coalescent_time = prev_coalescent_time + u;

        while ( next_coalescent_time < branch_length && j > 1 )
        {
            // randomly coalesce two lineages
            size_t index = static_cast<size_t>( floor(rng->uniform01()*initial_individuals_at_branch.size()) );
            TopologyNode *left = initial_individuals_at_branch[index];
            initial_individuals_at_branch.erase( initial_individuals_at_branch.begin() + index);

            index = static_cast<size_t>( floor(rng->uniform01()*initial_individuals_at_branch.size()) );
            TopologyNode *right = initial_individuals_at_branch[index];
            initial_individuals_at_branch.erase( initial_individuals_at_branch.begin() + index);

            TopologyNode *new_parent = new TopologyNode();
            new_parent->addChild(left);
            left->setParent(new_parent);
            new_parent->addChild(right);
            right->setParent(new_parent);

            root = new_parent;

            if ( root == NULL )
            {
                std::cerr << "Oh, the root is NULL :(" << std::endl;
            }

            initial_individuals_at_branch.push_back( new_parent );

            nodes_2_ages[new_parent] = next_coalescent_time + sp_node->getAge();


            prev_coalescent_time = next_coalescent_time;
            j--;
            n_pairs = j * (j-1) / 2.0;
            lambda = n_pairs * theta ;
            u = RbStatistics::Exponential::rv( lambda, *rng);
            next_coalescent_time = prev_coalescent_time + u;
        }

        if ( sp_parent_node != NULL )
        {
            std::vector<TopologyNode *> &incoming_lineages = individuals_per_branch[sp_parent_node];
            incoming_lineages.insert(incoming_lineages.end(), initial_individuals_at_branch.begin(), initial_individuals_at_branch.end());
        }


    }

    // the time tree object (topology + times)
    Tree *psi = new Tree();

    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );

    // initialize the topology by setting the root
    psi->setRoot(root, true);

    for ( std::map<TopologyNode*, double>::iterator it = nodes_2_ages.begin(); it != nodes_2_ages.end(); ++it)
    {
        TopologyNode *node = it->first;
        node->setAge( it->second );
    }

    // finally store the new value
    value = psi;

    resetTipAllocations();
}


/** Swap a parameter of the distribution */
void AbstractMultispeciesCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == species_tree )
    {
        species_tree = static_cast<const TypedDagNode< Tree >* >( newP );
    }

}
