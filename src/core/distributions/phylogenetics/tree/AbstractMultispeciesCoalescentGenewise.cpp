#include <stddef.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "AbstractMultispeciesCoalescentGenewise.h"
#include "DistributionExponential.h"
#include "ModelVector.h"
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

#include "RbMathFunctions.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

AbstractMultispeciesCoalescentGenewise::AbstractMultispeciesCoalescentGenewise(const TypedDagNode<Tree> *sp, RevBayesCore::RbVector< RevBayesCore::RbVector<Taxon> > t, size_t ngt) : TypedDistribution< RbVector<Tree> >( new RbVector<Tree>() ),
    taxa( t ),
    species_tree( sp ),
    num_taxa( ),
    log_tree_topology_prob ( 0.0 ),
    num_gene_trees( ngt ),
    gene_trees( )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( species_tree );

    // Get num_taxa for all gene trees
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        num_taxa.push_back( taxa[i].size() );
    }

    // Get species names
    std::vector< std::set<std::string> > species_names;
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        std::set<std::string> sn;

        for (RevBayesCore::RbIterator<RevBayesCore::Taxon> it=taxa[i].begin(); it!=taxa[i].end(); ++it)
        {
            sn.insert( it->getSpeciesName() );
        }

        species_names.push_back( sn );
    }

    // Get combinatorial topology prob
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        double ln_fact = RbMath::lnFactorial((int)(num_taxa[i]));
        log_tree_topology_prob += (num_taxa[i] - 1) * RbConstants::LN2 - 2.0 * ln_fact - std::log( num_taxa[i] ) ;
    }

    // Clear the current values
    //value->clear();
    //individuals_per_branch_genewise.clear();

    //std::cout << "num gene trees: " << gene_trees.size() << std::endl;

    gene_trees.clear();

    redrawValue();

}


AbstractMultispeciesCoalescentGenewise::~AbstractMultispeciesCoalescentGenewise()
{

}


void AbstractMultispeciesCoalescentGenewise::attachTimes(std::vector<Tree*> psi, std::vector< std::vector<TopologyNode *> > &tips, size_t index, const std::vector< std::vector<double> > &times) {

    for (size_t i=0; i<num_gene_trees; ++i)
    {
        if (index < num_taxa[i]-1)
        {
            // Get the rng
            RandomNumberGenerator* rng = GLOBAL_RNG;

            // randomly draw one node from the list of tips
            size_t tip_index = static_cast<size_t>( floor(rng->uniform01()*tips[i].size()) );

            // get the node from the list
            TopologyNode* parent = tips[i].at(tip_index);
            psi[i]->getNode( parent->getIndex() ).setAge( times[i][num_taxa[i] - index - 2] );

            // remove the randomly drawn node from the list
            tips[i].erase(tips[i].begin()+tip_index);

            // add a left child
            TopologyNode* leftChild = &parent->getChild(0);
            if ( !leftChild->isTip() )
            {
                tips[i].push_back(leftChild);
            }

            // add a right child
            TopologyNode* rightChild = &parent->getChild(1);
            if ( !rightChild->isTip() )
            {
                tips[i].push_back(rightChild);
            }

            // recursive call to this function
            attachTimes(psi, tips, index+1, times);
        }
    }
}


void AbstractMultispeciesCoalescentGenewise::buildRandomBinaryTree(std::vector< std::vector<TopologyNode*> > &tips)
{

    for (size_t i=0; i<num_gene_trees; ++i)
    {
        if (tips[i].size() < num_taxa[i])
        {
            // Get the rng
            RandomNumberGenerator* rng = GLOBAL_RNG;

            // randomly draw one node from the list of tips
            size_t index = static_cast<size_t>( floor(rng->uniform01()*tips[i].size()) );

            // get the node from the list
            TopologyNode* parent = tips[i].at(index);

            // remove the randomly drawn node from the list
            tips[i].erase(tips[i].begin()+index);

            // add a left child
            TopologyNode* leftChild = new TopologyNode(0);
            parent->addChild(leftChild);
            leftChild->setParent(parent);
            tips[i].push_back(leftChild);

            // add a right child
            TopologyNode* rightChild = new TopologyNode(0);
            parent->addChild(rightChild);
            rightChild->setParent(parent);
            tips[i].push_back(rightChild);

            // recursive call to this function
            buildRandomBinaryTree(tips);
        }
    }

}



double AbstractMultispeciesCoalescentGenewise::computeLnProbability( void )
{
    resetTipAllocations();

    // variable declarations and initialization
    double ln_prob_coal = 0;

    const Tree &sp = species_tree->getValue();

    ln_prob_coal = recursivelyComputeLnProbability( sp.getRoot() );

    // std::cout << sp.getNewickRepresentation() << std::endl;
    // std::cout << "final ln prob coal: " << ln_prob_coal << std::endl;

    return ln_prob_coal; // + logTreeTopologyProb;

}



double AbstractMultispeciesCoalescentGenewise::drawNe( size_t index )
{
    return 1.0;
}



double AbstractMultispeciesCoalescentGenewise::recursivelyComputeLnProbability( const RevBayesCore::TopologyNode &species_node )
{
    double ln_prob_coal = 0;

    if ( species_node.isTip() == false )
    {
        for (size_t i=0; i<num_gene_trees; i++)
        {
            individuals_per_branch_genewise[ i ][ species_node.getIndex() ].clear();
        }

        for (size_t j=0; j<species_node.getNumberOfChildren(); ++j)
        {
            ln_prob_coal += recursivelyComputeLnProbability( species_node.getChild(j) );
        }
    }

    double species_age = species_node.getAge();
    double parent_species_age = RbConstants::Double::inf;

    if ( species_node.isRoot() == false )
    {
        const TopologyNode &species_parent_node = species_node.getParent();
        parent_species_age = species_parent_node.getAge();
    }

    std::vector< std::vector<double> > coal_times_genewise;
    std::vector<size_t> initial_individuals_sizes_genewise;
    std::vector< std::set<const TopologyNode*> > remaining_individuals_genewise;

    for (size_t i=0; i<num_gene_trees; i++)
    {
        // create a local copy of the individuals per branch
        const std::set<const TopologyNode*> &current_initial_individuals = individuals_per_branch_genewise[i][species_node.getIndex()];
        std::set<const TopologyNode*> current_remaining_individuals = current_initial_individuals;

        // Get number of initial individuals for this gene in this branch
        initial_individuals_sizes_genewise.push_back( current_initial_individuals.size() );

        // Get all coalescent events among the individuals for this gene
        std::vector<double> current_coal_times;
        std::map<double, const TopologyNode *> current_coal_times_2_nodes;

        for ( std::set<const TopologyNode*>::iterator it = current_remaining_individuals.begin(); it != current_remaining_individuals.end(); ++it)
        {
            const TopologyNode *ind = (*it);
            if ( ind->isRoot() == false )
            {
                const TopologyNode &parent = ind->getParent();
                double parent_age = parent.getAge();
                current_coal_times_2_nodes[ parent_age ] = &parent;
            }
        }

        double current_time = species_age;

        while ( current_time < parent_species_age && current_coal_times_2_nodes.size() > 0 )
        {
            const TopologyNode *parent = current_coal_times_2_nodes.begin()->second;
            double parent_age = parent->getAge();
            current_time = parent_age;

            if ( parent_age < parent_species_age )
            { //Coalescence in the species tree branch

                // get the left and right child of the parent
                const TopologyNode *left = &parent->getChild( 0 );
                const TopologyNode *right = &parent->getChild( 1 );
                if ( current_remaining_individuals.find( left ) == current_remaining_individuals.end() || current_remaining_individuals.find( right ) == current_remaining_individuals.end() )
                {
                    // one of the children does not belong to this species tree branch
                    return RbConstants::Double::neginf;
                }

                //We remove the coalescent event and the coalesced lineages
                current_coal_times_2_nodes.erase( current_coal_times_2_nodes.begin() );
                current_remaining_individuals.erase( current_remaining_individuals.find( left ) );
                current_remaining_individuals.erase( current_remaining_individuals.find( right ) );

                //We insert the parent in the vector of lineages in this branch
                current_remaining_individuals.insert( parent );
                if ( parent->isRoot() == false )
                {
                    const TopologyNode *grand_parent = &parent->getParent();
                    current_coal_times_2_nodes[ grand_parent->getAge() ] = grand_parent;
                }

                current_coal_times.push_back( parent_age );


            } //End if coalescence in the species tree branch
            else
            { //No more coalescences in this species tree branch

                // jump out of the while loop
                //                currentTime = speciesAge;
                break;
            }

        } // end of while loop

        // Add coal times and remaining individuals to genewise data structures
        coal_times_genewise.push_back( current_coal_times );
        remaining_individuals_genewise.push_back( current_remaining_individuals );

        // Merge the two sets of individuals that go into the next species
        if ( species_node.isRoot() == false )
        {
            std::set<const TopologyNode *> &current_incoming_lineages = individuals_per_branch_genewise[ i ][ species_node.getParent().getIndex() ];
            current_incoming_lineages.insert( remaining_individuals_genewise[i].begin(), remaining_individuals_genewise[i].end());
        }
    }

    // Calculate log likelihood for the branch
    ln_prob_coal += computeLnCoalescentProbability(initial_individuals_sizes_genewise, coal_times_genewise, species_age, parent_species_age, species_node.getIndex(), species_node.isRoot() == false);

    // // Merge the two sets of individuals that go into the next species
    // for (size_t i=0; i<num_gene_trees; ++i)
    // {
    //     if ( species_node.isRoot() == false )
    //     {
    //         std::set<const TopologyNode *> &current_incoming_lineages = individuals_per_branch_genewise[ i ][ species_node.getParent().getIndex() ];
    //         current_incoming_lineages.insert( remaining_individuals_genewise[i].begin(), remaining_individuals_genewise[i].end());
    //     }
    // }

    return ln_prob_coal;
}



void AbstractMultispeciesCoalescentGenewise::redrawValue( void )
{

    simulateTrees();

}


/**
 * Build the map between species nodes and embedded individuals for each gene
 * (i.e., individuals_per_branch_genewise).
 */
void AbstractMultispeciesCoalescentGenewise::resetTipAllocations( void )
{
    individuals_per_branch_genewise.clear();

    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();

    // First let's create a map from species names to the tip nodes
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            const std::string &name = (*it)->getName();
            species_names_2_species_nodes[name] = *it;
        }
    }

    // Now we build the map between the species nodes and the embedded individuals

    // Get reference to gene trees
    std::vector<Tree>& values = *value;

    for (size_t i=0; i<num_gene_trees; ++i)
    {
        Tree current_tree = values[i];

        // Create a map between individual names and species names
        std::map<std::string, std::string> individual_names_2_species_names;
        for (RevBayesCore::RbIterator<RevBayesCore::Taxon> it=taxa[i].begin(); it!=taxa[i].end(); ++it)
        {
            const std::string &name = it->getName();
            individual_names_2_species_names[name] = it->getSpeciesName();
        }

        // Create maps for the individuals to branches
        individuals_per_branch_genewise.push_back( std::vector< std::set< const TopologyNode* > >(sp.getNumberOfNodes(), std::set< const TopologyNode* >() ) );

        for (size_t j=0; j<num_taxa[i]; ++j)
        {
            const TopologyNode &n = values[i].getNode( j );
            const std::string &individual_name = n.getName();
            const std::string &species_name = individual_names_2_species_names[ individual_name ];

            TopologyNode *species_node = species_names_2_species_nodes[species_name];
            individuals_per_branch_genewise[ i ][ species_node->getIndex() ].insert( &n );
        }
    }
}



/**
 * Set the current value.
 */
void AbstractMultispeciesCoalescentGenewise::setValue(RbVector<Tree>* v, bool f )
{

    // free memory
    if (value != v)
    {
        delete value;
    }

    value = v;

    resetTipAllocations();

}



void AbstractMultispeciesCoalescentGenewise::simulateTrees( void )
{
    // Clear the current values
    value->clear();
    individuals_per_branch_genewise.clear();

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Get the species tree and create a map from species names to the nodes of the species tree
    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    std::map<std::string, TopologyNode * > species_names_2_nodes;

    // std::cout << "sp tree: " << sp << std::endl;

    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            // std::cout << "name: " << name << std::endl;
            species_names_2_nodes[name] = *it;
        }
    }

    // Simulate all gene trees
    // Set up genewise data structures
    std::vector< std::map< const TopologyNode *, std::vector< TopologyNode* > > > individuals_per_branch_genewise;

    // We loop through each gene tree to get the information about individuals
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        std::map< const TopologyNode *, std::vector< TopologyNode* > > current_individuals_per_branch;

        // First we deal with the tips
        for (RevBayesCore::RbIterator<RevBayesCore::Taxon> it=taxa[i].begin(); it!=taxa[i].end(); ++it)
        {
            TopologyNode *n = new TopologyNode( *it );
            const std::string &species_name = n->getSpeciesName();

            if ( species_name == "" )
            {
                throw RbException("Cannot match a taxon without species to a tip in the species tree. The taxon map is probably wrong.");
            }

            TopologyNode *species_node = species_names_2_nodes[species_name];

            if ( species_node == NULL )
            {
                throw RbException("Could not match a taxon with name" + species_name + " to any of the tips in the species tree.");
            }

            n->setAge( 0.0 ); // Assumes all taxa are modern

            // Add nodes for this branch
            std::vector<TopologyNode * > &current_nodes_at_node = current_individuals_per_branch[ species_node ];
            current_nodes_at_node.push_back( n );
        }

        // // Keep track of individuals per branch
        // individuals_per_branch_genewise.push_back(current_individuals_per_branch);

        // Now deal with the internal nodes
        std::map<TopologyNode *, double> current_nodes_2_ages;
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

            std::vector<TopologyNode * > current_initial_individuals_at_branch = current_individuals_per_branch[sp_node];
            double branch_ne = drawNe( sp_node->getIndex() );

            double theta = 2.0 / branch_ne;

            double prev_coalescent_time = 0.0;

            size_t j = current_initial_individuals_at_branch.size();
            double n_pairs = j * (j-1) / 2.0;
            double lambda = n_pairs * theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);
            double next_coalescent_time = prev_coalescent_time + u;

            while ( next_coalescent_time < branch_length && j > 1 )
            {
                // randomly coalesce two lineages
                size_t index = static_cast<size_t>( floor(rng->uniform01()*current_initial_individuals_at_branch.size()) );
                TopologyNode *left = current_initial_individuals_at_branch[index];
                current_initial_individuals_at_branch.erase( current_initial_individuals_at_branch.begin() + index);

                index = static_cast<size_t>( floor(rng->uniform01()*current_initial_individuals_at_branch.size()) );
                TopologyNode *right = current_initial_individuals_at_branch[index];
                current_initial_individuals_at_branch.erase( current_initial_individuals_at_branch.begin() + index);

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

                current_initial_individuals_at_branch.push_back( new_parent );

                current_nodes_2_ages[new_parent] = next_coalescent_time + sp_node->getAge();


                prev_coalescent_time = next_coalescent_time;
                j--;
                n_pairs = j * (j-1) / 2.0;
                lambda = n_pairs * theta ;
                u = RbStatistics::Exponential::rv( lambda, *rng);
                next_coalescent_time = prev_coalescent_time + u;
            }

            if ( sp_parent_node != NULL )
            {
                std::vector<TopologyNode *> &current_incoming_lineages = current_individuals_per_branch[sp_parent_node];
                current_incoming_lineages.insert(current_incoming_lineages.end(), current_initial_individuals_at_branch.begin(), current_initial_individuals_at_branch.end());
            }
        }

        // Now save the individuals per branch for this gene
        individuals_per_branch_genewise.push_back(current_individuals_per_branch);

        // the time tree object (topology + times)
        Tree *psi = new Tree();

        // internally we treat unrooted topologies the same as rooted
        psi->setRooted( true );

        // initialize the topology by setting the root
        psi->setRoot(root, true);

        for ( std::map<TopologyNode*, double>::iterator it = current_nodes_2_ages.begin(); it != current_nodes_2_ages.end(); ++it)
        {
            TopologyNode *node = it->first;
            node->setAge( it->second );
        }

        // finally store the new value
        value->push_back(*psi);
        gene_trees.push_back(psi);

    }

    resetTipAllocations();

}



std::vector<Tree*> AbstractMultispeciesCoalescentGenewise::getTrees(void) const
{
    return gene_trees;
}



size_t AbstractMultispeciesCoalescentGenewise::getNumberOfGeneTrees(void) const
{
    return num_gene_trees;
}



/** Swap a parameter of the distribution */
void AbstractMultispeciesCoalescentGenewise::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == species_tree )
    {
        species_tree = static_cast<const TypedDagNode< Tree >* >( newP );
    }

    // if ( oldP == gene_trees )
    // {
    //     gene_trees = static_cast< std::vector< Tree* > >( newP );
    // }

}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const AbstractMultispeciesCoalescentGenewise& x)
{

    std::vector<Tree*> trees = x.getTrees();

    for ( size_t i=0; i<x.getNumberOfGeneTrees(); ++i )
    {
        o << trees[i]->getNewickRepresentation();
    }

    return o;
}
