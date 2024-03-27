#include <cstddef>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

#include "DistributionBeta.h"
#include "SpeciesSubtreeScaleBetaProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TreeUtilities.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SpeciesSubtreeScaleBetaProposal::SpeciesSubtreeScaleBetaProposal( StochasticNode<Tree> *sp, double a ) : Proposal(),
    speciesTree( sp ),
    geneTrees( ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( speciesTree );

}


/**
 * Add a new DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesSubtreeScaleBetaProposal::addGeneTree(StochasticNode<Tree> *gt)
{
    // check if this node isn't already in our list
    bool exists = false;
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            exists = true;
            break;
        }
    }

    // only add this variable if it doesn't exist in our list already
    if ( exists == false )
    {
        geneTrees.push_back( gt );
        addNode( gt );
    }

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SpeciesSubtreeScaleBetaProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SpeciesSubtreeScaleBetaProposal* SpeciesSubtreeScaleBetaProposal::clone( void ) const
{

    return new SpeciesSubtreeScaleBetaProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SpeciesSubtreeScaleBetaProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesSubtreeScaleBeta";

    return name;
}


double SpeciesSubtreeScaleBetaProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double SpeciesSubtreeScaleBetaProposal::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = speciesTree->getValue();

    // pick a random node which is not the root and not a tip
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );

    TopologyNode& parent = node->getParent();

    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();

    // now we store all necessary values
    storedNode = node;
    storedSpeciesTreeAges = std::vector<double>(tau.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(*node, storedSpeciesTreeAges);

    // lower bound
    double min_age = TreeUtilities::getOldestTipAge(*node);

    // draw new ages
    double current_value = my_age / (parent_age - min_age);
    double a = alpha + 1.0;
    double b = (a-1.0) / current_value - a + 2.0;
    double new_value = RbStatistics::Beta::rv(a, b, *rng);

    // Sebastian: This is for debugging to test if the proposal's acceptance rate is 1.0 as it should be!
//    new_value = current_value;

    double my_new_age = new_value * (parent_age - min_age) + min_age;

    double scaling_factor = my_new_age / my_age;

    size_t num_nodes = node->getNumberOfNodesInSubtree( false );

    storedGeneTreeAges = std::vector<std::vector<double> >(geneTrees.size(), std::vector<double>());

    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& gene_tree = geneTrees[i]->getValue();

        std::vector<TopologyNode*> nodes = getOldestNodesInPopulation(gene_tree, *node );

        storedGeneTreeAges[i] = std::vector<double>(gene_tree.getNumberOfNodes(), 0.0);
        TreeUtilities::getAges(gene_tree.getRoot(), storedGeneTreeAges[i]);

        for (size_t j=0; j<nodes.size(); ++j)
        {
            // add the number of nodes that we are going to scale in the subtree
            num_nodes += nodes[j]->getNumberOfNodesInSubtree( false );

            // rescale the subtree of this gene tree
            TreeUtilities::rescaleSubtree(*nodes[j], scaling_factor );

        }

    }

    // Sebastian: We need to work on a mechanism to make these proposal safe for non-ultrametric trees!
    //    if (min_age != 0.0)
    //    {
    //        for (size_t i = 0; i < tau.getNumberOfTips(); i++)
    //        {
    //            if (tau.getNode(i).getAge() < 0.0)
    //            {
    //                return RbConstants::Double::neginf;
    //            }
    //        }
    //    }

    // rescale the subtree of the species tree
    TreeUtilities::rescaleSubtree(*node, scaling_factor );

    // compute the Hastings ratio
    double forward = RbStatistics::Beta::lnPdf(a, b, new_value);
    double new_a = alpha + 1.0;
    double new_b = (new_a-1.0) / new_value - new_a + 2.0;
    double backward = RbStatistics::Beta::lnPdf(new_a, new_b, current_value);

    double lnHastingsratio = (num_nodes > 1 ? ( log(scaling_factor) + backward - forward ) * (num_nodes-1) : 0.0 );

    return lnHastingsratio;
}


std::vector<TopologyNode*> SpeciesSubtreeScaleBetaProposal::getOldestNodesInPopulation( Tree &tau, TopologyNode &n )
{

    // I need all the oldest nodes/subtrees that have the same tips.
    // Those nodes need to be scaled too.

    // get the beginning and ending age of the population
    double max_age = -1.0;
    if ( n.isRoot() == false )
    {
        max_age = n.getParent().getAge();
    }

    // get all the taxa from the species tree that are descendants of node i
    std::vector<TopologyNode*> speciesTaxa;
    TreeUtilities::getTaxaInSubtree( n, speciesTaxa );

    // get all the individuals
    std::set<TopologyNode*> individualTaxa;
    for (size_t i = 0; i < speciesTaxa.size(); ++i)
    {
        const std::string &name = speciesTaxa[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individualTaxa.insert( ind[j] );
        }
    }

    // create the set of the nodes within this population
    std::set<TopologyNode*> nodesInPopulationSet;

    // now go through all nodes in the gene
    while ( individualTaxa.empty() == false )
    {
        // get the first element
        std::set<TopologyNode*>::iterator it = individualTaxa.begin();

        // store the pointer
        TopologyNode *geneNode = *it;

        // and now remove the element from the list
        individualTaxa.erase( it );

        // add this node to our list of node we need to scale, if:
        // a) this is the root node
        // b) this is not the root and the age of the parent node is larger than the parent's age of the species node
        if ( geneNode->isRoot() == false && ( max_age == -1.0 || max_age > geneNode->getParent().getAge() ) )
        {
            // push the parent to our current list
            individualTaxa.insert( &geneNode->getParent() );
        }
        else if ( geneNode->isTip() == false )
        {
            // add this node if it is within the age of our population
            nodesInPopulationSet.insert( geneNode );
        }

    }

    // convert the set into a vector
    std::vector<TopologyNode*> nodesInPopulation;
    for (std::set<TopologyNode*>::iterator it = nodesInPopulationSet.begin(); it != nodesInPopulationSet.end(); ++it)
    {
        nodesInPopulation.push_back( *it );
    }

    return nodesInPopulation;
}


/**
 *
 */
void SpeciesSubtreeScaleBetaProposal::prepareProposal( void )
{

}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SpeciesSubtreeScaleBetaProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }

}


/**
 * Remove a DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesSubtreeScaleBetaProposal::removeGeneTree(StochasticNode<Tree> *gt)
{
    // remove it from our list
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            geneTrees.erase( geneTrees.begin() + i );
            --i;
        }
    }

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SpeciesSubtreeScaleBetaProposal::undoProposal( void )
{
    // undo the proposal
    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& geneTree = geneTrees[i]->getValue();

        TreeUtilities::setAges(geneTree.getRoot(), storedGeneTreeAges[i]);

    }
    TreeUtilities::setAges(*storedNode, storedSpeciesTreeAges);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SpeciesSubtreeScaleBetaProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if ( oldN == speciesTree )
    {
        speciesTree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else
    {
        for ( size_t i=0; i<geneTrees.size(); ++i )
        {
            if ( oldN == geneTrees[i] )
            {
                geneTrees[i] = static_cast<StochasticNode<Tree>* >(newN) ;
            }
        }
    }

}


void SpeciesSubtreeScaleBetaProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SpeciesSubtreeScaleBetaProposal::tune( double rate )
{

    if ( rate > 0.234 )
    {
        alpha /= (1.0 + ((rate-0.234)/0.766) );
    }
    else
    {
        alpha *= (2.0 - rate/0.234 );
    }

}
