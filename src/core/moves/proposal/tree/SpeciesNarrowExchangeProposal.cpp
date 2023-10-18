#include "SpeciesNarrowExchangeProposal.h"

#include <cstddef>
#include <cmath>
#include <iterator>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
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
SpeciesNarrowExchangeProposal::SpeciesNarrowExchangeProposal( StochasticNode<Tree> *sp ) : Proposal(),
    speciesTree( sp ),
    geneTrees( )
{
    // tell the base class to add the node
    addNode( speciesTree );

}


/**
 * Add a new DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesNarrowExchangeProposal::addGeneTree(StochasticNode<Tree> *gt)
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
void SpeciesNarrowExchangeProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SpeciesNarrowExchangeProposal* SpeciesNarrowExchangeProposal::clone( void ) const
{

    return new SpeciesNarrowExchangeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SpeciesNarrowExchangeProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesNarrowExchange";

    return name;
}


double SpeciesNarrowExchangeProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double SpeciesNarrowExchangeProposal::doProposal( void )
{

    // empty the previous vectors
    storedGeneTreeNodes.clear();
    storedOldBrothers.clear();

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = speciesTree->getValue();

    // pick a random node which is not the root and neithor a direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );

    TopologyNode& parent = node->getParent();
    TopologyNode& grandparent = parent.getParent();
    TopologyNode* uncle = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == &parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    TopologyNode* brother = &parent.getChild( 0 );
    // check if we got the correct child
    if ( brother == node )
    {
        brother = &parent.getChild( 1 );
    }

    // we need to work with the times
    double parent_age   = parent.getAge();
    double uncles_age   = uncle->getAge();

    if ( uncles_age < parent_age )
    {
        failed = false;

        double ln_hastings_ratio = 0.0;

        // now we store all necessary values
        storedChoosenNode   = node;
        storedUncle         = uncle;

        // now we need to find for each gene tree the nodes that need to be moved as well
        // only nodes that have a coalescent event within the lifetime of the parents populations
        // from lineages belonging to the chosen node with lineages belonging to the brother population
        // need to be changed
        for ( size_t i=0; i<geneTrees.size(); ++i )
        {
            // get the i-th gene tree
            Tree& geneTree = geneTrees[i]->getValue();

            std::vector<TopologyNode*> nodes = getNodesToChange(geneTree, *node, *brother );

            // get the set of nodes in my uncles populations
            // these are the nodes that are possible re-attachment points
            std::set<TopologyNode*> new_siblings = getOldestSubtreesNodesInPopulation(geneTree, *uncle);
            std::set<TopologyNode*> old_siblings = getOldestSubtreesNodesInPopulation(geneTree, *brother);

            for (size_t j=0; j<nodes.size(); ++j)
            {

                TopologyNode *the_gene_node = nodes[i];

                // first we need to compute the backward probability
                std::set<TopologyNode*> old_candidate_siblings = getPossibleSiblings(the_gene_node, old_siblings);

                // add the backward probability to the hastings ratio
                ln_hastings_ratio += log( old_siblings.size() );

                // then we need to compute the forward probability
                std::set<TopologyNode*> new_candidate_siblings = getPossibleSiblings(the_gene_node, new_siblings);

                // add the forward  probability to the hastings ratio
                ln_hastings_ratio += log( new_candidate_siblings.size() );

                // actually pick a new sibling
                size_t new_index = size_t( floor(rng->uniform01() * new_candidate_siblings.size() ) );
                std::set<TopologyNode*>::iterator it = new_candidate_siblings.begin();
                std::advance(it,new_index);
                TopologyNode *new_child = *it;

                // store nodes
                storedGeneTreeNodes.push_back( the_gene_node );
                TopologyNode &the_parent = the_gene_node->getParent();
                TopologyNode *old_brother = &the_parent.getChild( 0 );
                if ( old_brother == the_gene_node )
                {
                    old_brother = &the_parent.getChild( 1 );
                }
                storedOldBrothers.push_back( old_brother );

                // perform a prune and regraft move
                prune( &the_parent, the_gene_node );
                regraft( the_gene_node, new_child );

            }

            // Sebastian: This is only for debugging. It makes the code slower. Hopefully it is not necessary anymore.
            //        geneTrees[i]->touch( true );

        }


        // now exchange the two nodes
        grandparent.removeChild( uncle );
        parent.removeChild( node );
        grandparent.addChild( node );
        parent.addChild( uncle );
        node->setParent( &grandparent );
        uncle->setParent( &parent );

        return 0.0;
//        return ln_hastings_ratio;
    }
    else
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

}


std::vector<TopologyNode*> SpeciesNarrowExchangeProposal::getNodesToChange( Tree &tau, TopologyNode &node, TopologyNode &brother )
{
    // get the age of the parent node
    double parent_age = node.getParent().getAge();

    // get the nodes representing the subtrees in my population
    std::set<TopologyNode*> my_nodes = getOldestSubtreesNodesInPopulation(tau, node);

    // get the nodes representing the subtrees in my population
    std::set<TopologyNode*> my_brothers_nodes = getNodesInPopulation(tau, brother);

    // create a vector for the nodes that need changing
    std::vector<TopologyNode*> nodes_to_change;

    // now we need to find the nodes that coalesce before the age of my parent
    while ( my_nodes.empty() == false )
    {
        // get the first node
        TopologyNode *the_node = *(my_nodes.begin());

        // remove the node from my list
        my_nodes.erase( my_nodes.begin() );

        // now get the parent of this node
        TopologyNode *the_parent = &the_node->getParent();

        // check the age
        if ( the_parent->getAge() < parent_age )
        {
            // now get the brother of the node
            TopologyNode *the_brother = &the_parent->getChild(0);

            // check if we got the wrong child
            if ( the_brother == the_node )
            {
                the_brother = &the_parent->getChild(1);
            }

            // now check in which population the brother lived
            if ( my_nodes.find( the_brother ) != my_nodes.end() )
            {
                // we both lived in the same population
                // -> remove the brother as well and add the parent
                my_nodes.erase( the_brother );
                my_nodes.insert( the_parent );
            }
            else if ( my_brothers_nodes.find( the_brother ) != my_brothers_nodes.end() )
            {
                // we found a node which we need to change
                nodes_to_change.push_back( the_node );
            }
            else
            {
                // this node doesn't coalesce with a node from the brothers population
            }

        }

    }

    return nodes_to_change;

}


std::set<TopologyNode*> SpeciesNarrowExchangeProposal::getNodesInPopulation( Tree &tau, TopologyNode &n )
{

    // I need all the oldest nodes/subtrees that have the same tips.
    // Those nodes need to be scaled too.

    // get the beginning and ending age of the population
    double max_age = -1.0;
    if ( n.isRoot() == false )
    {
        max_age = n.getParent().getAge();
    }

    // get all the taxa from the species tree that are descendants of the node
    double min_age = n.getAge();
    std::vector<TopologyNode*> species_taxa;
    TreeUtilities::getTaxaInSubtree( n, species_taxa );

    // get all the individuals
    std::set<TopologyNode*> individual_taxa;
    for (size_t i = 0; i < species_taxa.size(); ++i)
    {
        const std::string &name = species_taxa[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individual_taxa.insert( ind[j] );
        }
    }

    // create the set of the nodes within this population
    std::set<TopologyNode*> nodesInPopulationSet;

    // now go through all nodes in the gene
    while ( individual_taxa.empty() == false )
    {
        std::set<TopologyNode*>::iterator it = individual_taxa.begin();
        individual_taxa.erase( it );

        TopologyNode *gene_node = *it;

        if ( gene_node->getAge() < min_age )
        {
            // the age of the node is younger than the populations start age
            // -> add the node to the current working set
            individual_taxa.insert( &gene_node->getParent() );

        }
        else if ( gene_node->getAge() < max_age || max_age == -1.0 )
        {
            // the age of the node is within the population
            // -> add the node to the current return set
            nodesInPopulationSet.insert( gene_node );

            // if this is not the root then we need to add the parent node to the working set
            if ( gene_node->isRoot() == false )
            {
                individual_taxa.insert( &gene_node->getParent() );

            }

        }

    }

    return nodesInPopulationSet;
}


std::set<TopologyNode*> SpeciesNarrowExchangeProposal::getOldestSubtreesNodesInPopulation( Tree &tau, TopologyNode &n )
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
    std::vector<TopologyNode*> species_taxa;
    TreeUtilities::getTaxaInSubtree( n, species_taxa );

    // get all the individuals
    std::set<TopologyNode*> individualTaxa;
    for (size_t i = 0; i < species_taxa.size(); ++i)
    {
        const std::string &name = species_taxa[i]->getName();
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
        std::set<TopologyNode*>::iterator it = individualTaxa.begin();
        individualTaxa.erase( it );

        TopologyNode *geneNode = *it;

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

    return nodesInPopulationSet;
}


std::set<TopologyNode*> SpeciesNarrowExchangeProposal::getPossibleSiblings(TopologyNode *n, const std::set<TopologyNode *> candidates)
{

    std::set<TopologyNode*> possible_new_siblings;

    for (std::set<TopologyNode*>::iterator it = candidates.begin(); it != candidates.end(); ++it)
    {
        TopologyNode *the_new_sibling = *it;
        TopologyNode &the_new_parent = the_new_sibling->getParent();

        double my_parents_age = n->getParent().getAge();

        if ( my_parents_age > the_new_sibling->getAge() && my_parents_age < the_new_parent.getAge() )
        {
            possible_new_siblings.insert( the_new_sibling );
        }

    }

    return possible_new_siblings;
}


/**
 *
 */
void SpeciesNarrowExchangeProposal::prepareProposal( void )
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
void SpeciesNarrowExchangeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    // no parameters

}



void SpeciesNarrowExchangeProposal::prune( TopologyNode *node, TopologyNode *child )
{
    // get the parent and brother node
    TopologyNode &parent = node->getParent();
    TopologyNode *brother = &node->getChild( 0 );
    if ( brother == child )
    {
        brother = &node->getChild( 1 );
    }
    // prune
    parent.removeChild( node );
    node->removeChild( brother );
    parent.addChild( brother );
    brother->setParent( &parent );

}



void SpeciesNarrowExchangeProposal::regraft( TopologyNode *node, TopologyNode *child )
{

    // regraft
    TopologyNode* newParent = &child->getParent();
    newParent->removeChild( child );
    newParent->addChild( node );
    node->setParent( newParent );
    node->addChild( child );
    child->setParent( node );
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SpeciesNarrowExchangeProposal::undoProposal( void )
{
    // we undo the proposal only if it didn't fail
    if ( !failed )
    {
        // undo the proposal
        TopologyNode& parent = storedUncle->getParent();
        TopologyNode& grandparent = storedChoosenNode->getParent();

        // now exchange the two nodes
        grandparent.removeChild( storedChoosenNode );
        parent.removeChild( storedUncle );
        grandparent.addChild( storedUncle );
        parent.addChild( storedChoosenNode );
        storedUncle->setParent( &grandparent );
        storedChoosenNode->setParent( &parent );

        // now also undo all the gene tree changes!
        for (size_t i = 0; i < storedGeneTreeNodes.size(); ++i)
        {
            TopologyNode *the_gene_node = storedGeneTreeNodes[i];
            TopologyNode *old_brother   = storedOldBrothers[i];

            // perform a prune and regraft move
            prune( &the_gene_node->getParent(), the_gene_node );
            regraft( the_gene_node, old_brother );
        }

    }


}


/**
 * Remove a DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesNarrowExchangeProposal::removeGeneTree(StochasticNode<Tree> *gt)
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
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void SpeciesNarrowExchangeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


void SpeciesNarrowExchangeProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SpeciesNarrowExchangeProposal::tune( double rate )
{

    // nothing to tune

}
