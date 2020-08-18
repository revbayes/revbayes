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

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

AbstractMultispeciesCoalescentGenewise::AbstractMultispeciesCoalescentGenewise(const TypedDagNode<Tree> *sp, RevBayesCore::RbVector< RevBayesCore::RbVector<Taxon> > t, size_t ngt) : TypedDistribution<Tree>( NULL ),
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

    // Get combinatorial topology probs
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        double ln_fact = RbMath::lnFactorial((int)(num_taxa[i]));
        log_tree_topology_prob += (num_taxa[i] - 1) * RbConstants::LN2 - 2.0 * ln_fact - std::log( num_taxa[i] ) ;
    }

    redrawValue();

}


AbstractMultispeciesCoalescentGenewise::~AbstractMultispeciesCoalescentGenewise()
{

}



// void AbstractMultispeciesCoalescentGenewise::attachTimes(std::vector<Tree *psi>, std::vector< std::vector<TopologyNode *> > &tips, size_t index, const std::vector< std::vector<double> > &times) {
//
//     for (size_t i=0; i<num_gene_trees; ++i)
//     {
//         if (index < num_taxa[i]-1)
//         {
//             // Get the rng
//             RandomNumberGenerator* rng = GLOBAL_RNG;
//
//             // randomly draw one node from the list of tips
//             size_t tip_index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
//
//             // get the node from the list
//             TopologyNode* parent = tips.at(tip_index);
//             psi->getNode( parent->getIndex() ).setAge( times[num_taxa - index - 2] );
//
//             // remove the randomly drawn node from the list
//             tips.erase(tips.begin()+tip_index);
//
//             // add a left child
//             TopologyNode* leftChild = &parent->getChild(0);
//             if ( !leftChild->isTip() )
//             {
//                 tips.push_back(leftChild);
//             }
//
//             // add a right child
//             TopologyNode* rightChild = &parent->getChild(1);
//             if ( !rightChild->isTip() )
//             {
//                 tips.push_back(rightChild);
//             }
//
//             // recursive call to this function
//             attachTimes(psi, tips, index+1, times);
//         }
//     }
// }
//
//
// void AbstractMultispeciesCoalescentGenewise::buildRandomBinaryTree(std::vector<TopologyNode*> &tips)
// {
//
//     if (tips.size() < num_taxa)
//     {
//         // Get the rng
//         RandomNumberGenerator* rng = GLOBAL_RNG;
//
//         // randomly draw one node from the list of tips
//         size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
//
//         // get the node from the list
//         TopologyNode* parent = tips.at(index);
//
//         // remove the randomly drawn node from the list
//         tips.erase(tips.begin()+index);
//
//         // add a left child
//         TopologyNode* leftChild = new TopologyNode(0);
//         parent->addChild(leftChild);
//         leftChild->setParent(parent);
//         tips.push_back(leftChild);
//
//         // add a right child
//         TopologyNode* rightChild = new TopologyNode(0);
//         parent->addChild(rightChild);
//         rightChild->setParent(parent);
//         tips.push_back(rightChild);
//
//         // recursive call to this function
//         buildRandomBinaryTree(tips);
//     }
// }


double AbstractMultispeciesCoalescentGenewise::computeLnProbability( void )
{
    resetTipAllocations();

    // variable declarations and initialization
    double ln_prob_coal = 0;

    const Tree &sp = species_tree->getValue();

    ln_prob_coal = recursivelyComputeLnProbability( sp.getRoot() );


    return ln_prob_coal; // + logTreeTopologyProb;

}


double AbstractMultispeciesCoalescentGenewise::drawNe( size_t index )
{

    return 1.0;
}


double AbstractMultispeciesCoalescentGenewise::recursivelyComputeLnProbability( const RevBayesCore::TopologyNode &species_node )
{

    double ln_prob_coal = 0;

    double species_age = species_node.getAge();
    double parent_species_age = RbConstants::Double::inf;

    std::vector< std::vector<double> > coal_times;
    std::vector< std::map<double, const TopologyNode *> > coal_times_2_nodes;
    std::vector< std::set<const TopologyNode*> > remaining_individuals;
    std::vector<size_t> initial_individuals_sizes;


    for (size_t i=0; i<num_gene_trees; ++i)
    {
        if ( species_node.isTip() == false )
        {
            individuals_per_branch_genewise[ i ][ species_node.getIndex() ].clear();

            for (size_t j=0; j<species_node.getNumberOfChildren(); ++j)
            {
                ln_prob_coal += recursivelyComputeLnProbability( species_node.getChild(j) );
            }
        }



        if ( species_node.isRoot() == false )
        {
            const TopologyNode &species_parent_node = species_node.getParent();
            parent_species_age = species_parent_node.getAge();
        }

        // create a local copy of the individuals per branch
        const std::set<const TopologyNode*> &initial_individuals = individuals_per_branch_genewise[i][species_node.getIndex()];
        remaining_individuals[i] = initial_individuals;

        // Get number of initial individuals
        initial_individuals_sizes[i] = initial_individuals.size();

        // get all coalescent events among the individuals
        for ( std::set<const TopologyNode*>::iterator it = remaining_individuals[i].begin(); it != remaining_individuals[i].end(); ++it)
        {
            const TopologyNode *ind = (*it);
            if ( ind->isRoot() == false )
            {
                const TopologyNode &parent = ind->getParent();
                double parent_age = parent.getAge();
                coal_times_2_nodes[ i ][ parent_age ] = &parent;
            }
        }

        double current_time = species_age;
        while ( current_time < parent_species_age && coal_times_2_nodes[i].size() > 0 )
        {

            const TopologyNode *parent = coal_times_2_nodes[i].begin()->second;
            double parent_age = parent->getAge();
            current_time = parent_age;

            if ( parent_age < parent_species_age )
            { //Coalescence in the species tree branch

                // get the left and right child of the parent
                const TopologyNode *left = &parent->getChild( 0 );
                const TopologyNode *right = &parent->getChild( 1 );
                if ( remaining_individuals[i].find( left ) == remaining_individuals[i].end() || remaining_individuals[i].find( right ) == remaining_individuals[i].end() )
                {
                    // one of the children does not belong to this species tree branch
                    return RbConstants::Double::neginf;
                }

                //We remove the coalescent event and the coalesced lineages
                coal_times_2_nodes[i].erase( coal_times_2_nodes[i].begin() );
                remaining_individuals[i].erase( remaining_individuals[i].find( left ) );
                remaining_individuals[i].erase( remaining_individuals[i].find( right ) );

                //We insert the parent in the vector of lineages in this branch
                remaining_individuals[i].insert( parent );
                if ( parent->isRoot() == false )
                {
                    const TopologyNode *grand_parent = &parent->getParent();
                    coal_times_2_nodes[ i ][ grand_parent->getAge() ] = grand_parent;
                }

                coal_times[i].push_back( parent_age );


            } //End if coalescence in the species tree branch
            else
            { //No more coalescences in this species tree branch

                // jump out of the while loop
                //                currentTime = speciesAge;
                break;
            }


        } // end of while loop
    }



    ln_prob_coal += computeLnCoalescentProbability(initial_individuals_sizes, coal_times, species_age, parent_species_age, species_node.getIndex(), species_node.isRoot() == false);



    // merge the two sets of individuals that go into the next species
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        if ( species_node.isRoot() == false )
        {
            std::vector< std::set<const TopologyNode *> > incoming_lineages;
            incoming_lineages[i] = individuals_per_branch_genewise[ i ][ species_node.getParent().getIndex() ];
            incoming_lineages[i].insert( remaining_individuals[i].begin(), remaining_individuals[i].end());
        }
    }


    return ln_prob_coal;
}


void AbstractMultispeciesCoalescentGenewise::redrawValue( void )
{

    simulateTrees();

}


void AbstractMultispeciesCoalescentGenewise::resetTipAllocations( void )
{

    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();

    // First let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            const std::string &name = (*it)->getName();
            species_names_2_species_nodes[name] = *it;
        }
    }

    // Second, let's create a map from individual names to the species names for each gene tree
    // This is a vector of maps between individual names and species names
    // std::vector< std::map<std::string, std::string> > individual_names_2_species_names_genewise;

    for (size_t i=0; i<num_gene_trees; ++i)
    {
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
            //const TopologyNode &n = value->getNode( j );
            const TopologyNode &n = gene_trees[i]->getNode( j );
            const std::string &individual_name = n.getName();
            const std::string &species_name = individual_names_2_species_names[ individual_name ];

            TopologyNode *species_node = species_names_2_species_nodes[species_name];
            individuals_per_branch_genewise[ i ][ species_node->getIndex() ].insert( &n );
        }
    }

}


// /**
//  * Set the current value.
//  */
// void AbstractMultispeciesCoalescentGenewise::setValues(std::vector<Tree*> *v, bool f )
// {
//
//     // delegate to super class
//     TypedDistribution<Tree>::setValue(v, f);
//
//     resetTipAllocations();
//
// }



void AbstractMultispeciesCoalescentGenewise::simulateTrees( void )
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Get the species tree and create a map from species names to the nodes of the species tree
    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    std::map<std::string, TopologyNode * > species_names_2_nodes;

    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            species_names_2_nodes[name] = *it;
        }
    }

    // Simulate all gene trees

    // Set up genewise data structures
    std::vector< std::map< const TopologyNode *, std::vector< TopologyNode* > > > individuals_per_branch_genewise;
    // std::vector< std::vector<TopologyNode * > > incoming_lineages_genewise;
    // std::vector< std::map<TopologyNode *, double> > nodes_2_ages_genewise;

    // We loop through each gene tree to get the information about individuals
    for (size_t i=0; i<num_gene_trees; ++i)
    {
        std::map< const TopologyNode *, std::vector< TopologyNode* > > current_individuals_per_branch;

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

            n->setAge( 0.0 );

            // Add nodes for this branch
            std::vector<TopologyNode * > &current_nodes_at_node = current_individuals_per_branch[ species_node ];
            current_nodes_at_node.push_back( n );
        }

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

            double theta = 1.0 / branch_ne;

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

        // Add individuals and lineages for this gene to the appropriate genewise vectors
        individuals_per_branch_genewise.push_back( current_individuals_per_branch );
        // incoming_lineages_genewise.push_back( current_incoming_lineages );
        // nodes_2_ages_genewise.push_back( current_nodes_2_ages );

        // finally store the new value
        //value[i] = *psi;
        gene_trees.push_back( psi );

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
