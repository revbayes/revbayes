#include "TreeNodeAgeUpdateProposal.h"

#include <cstddef>
#include <cmath>
#include <set>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
TreeNodeAgeUpdateProposal::TreeNodeAgeUpdateProposal( StochasticNode<Tree> *sp ) : Proposal(),
    speciesTree( sp ),
    geneTrees(  )
{
    // tell the base class to add the node
    addNode( speciesTree );

    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        addNode( geneTrees[i] );
    }

}


/**
 * Add a new DAG node holding a gene tree on which this move operates on.
 *
 */
void TreeNodeAgeUpdateProposal::addGeneTree(StochasticNode<Tree> *gt)
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
void TreeNodeAgeUpdateProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
TreeNodeAgeUpdateProposal* TreeNodeAgeUpdateProposal::clone( void ) const
{

    return new TreeNodeAgeUpdateProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& TreeNodeAgeUpdateProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesNodeTimeSlideUniform";

    return name;
}


double TreeNodeAgeUpdateProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double TreeNodeAgeUpdateProposal::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = speciesTree->getValue();
    
    // check that there is at least one node which is not the root, a tip, or the parent of a SA
    std::vector<TopologyNode*> eligible_nodes;
    
    for (auto& to_check: tau.getNodes())
    {
        if ( !to_check->isRoot() && !to_check->isTip() && !to_check->isSampledAncestorParent() )
        {
            eligible_nodes.push_back( to_check );
        }
    }
    
    if (eligible_nodes.size() == 0)
    {
        std::cerr << "mvSpeciesNodeTimeSlideUniform has no effect; the tree only contains the root, tips, and sampled ancestors." << std::endl;
        return RbConstants::Double::neginf;
    }

    // pick a random node which is not the root, a tip, or the parent of a SA
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() || node->isSampledAncestorParent() );

    TopologyNode& parent = node->getParent();

    // we need to work with the times: find the boundaries for the new age
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }

    // now we store all necessary values
    storedNode = node;
    storedAge = my_age;

    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;

    // Sebastian: This is for debugging to test if the proposal's acceptance rate is 1.0 as it should be!
//    my_new_age = my_age;

    int upslideNodes = 0;
    int downslideNodes = 0;

    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& geneTree = geneTrees[i]->getValue();

        std::vector<TopologyNode*> nodes = getNodesInPopulation(geneTree, *node );

        for (size_t j=0; j<nodes.size(); ++j)
        {

            double a = nodes[j]->getAge();
            double new_a = a;
            if ( a > my_age ) // gene coalescence was above speciation
            {
                ++upslideNodes;
                new_a = parent_age - (parent_age - my_new_age)/(parent_age - my_age) * (parent_age - a);
            }
            else
            {
              //std::cout << "##########################a < my_age"<<std::endl;
                ++downslideNodes;
                new_a = child_Age + (my_new_age - child_Age)/(my_age - child_Age) * (a - child_Age);
            }

            // set the new age of this gene tree node
            geneTree.getNode( nodes[j]->getIndex() ).setAge( new_a );
        }

        // Sebastian: This is only for debugging. It makes the code slower. Hopefully it is not necessary anymore.
//        geneTrees[i]->touch( true );

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


    // set the age of the species tree node
    tau.getNode( node->getIndex() ).setAge( my_new_age );

    // compute the Hastings ratio
    double lnHastingsratio = upslideNodes * log( (parent_age - my_new_age)/(parent_age - my_age) ) + downslideNodes * log( (my_new_age - child_Age)/(my_age - child_Age) );

  //  std::cout << "RUBBER BAND: parent_age: "<< parent_age << " my_age: "<< my_age << " my_new_age: "<< my_new_age << " child_Age: " << child_Age << " lnHastingsratio: " << lnHastingsratio<<std::endl;

    return lnHastingsratio;

}


std::vector<TopologyNode*> TreeNodeAgeUpdateProposal::getNodesInPopulation( Tree &tau, TopologyNode &n )
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
    double min_age_left = n.getChild(0).getAge();
    std::vector<TopologyNode*> speciesTaxa_left;
    TreeUtilities::getTaxaInSubtree( n.getChild(0), speciesTaxa_left );

    // get all the individuals
    std::set<TopologyNode*> individualTaxa_left;
    for (size_t i = 0; i < speciesTaxa_left.size(); ++i)
    {
        const std::string &name = speciesTaxa_left[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individualTaxa_left.insert( ind[j] );
        }
    }

    // create the set of the nodes within this population
    std::set<TopologyNode*> nodesInPopulationSet;

    // now go through all nodes in the gene
    while ( individualTaxa_left.empty() == false )
    {
        // get the first element
        std::set<TopologyNode*>::iterator it = individualTaxa_left.begin();

        // store the pointer
        TopologyNode *geneNode = *it;

        // and now remove the element from the list
        individualTaxa_left.erase( it );

        // add this node to our list of node we need to scale, if:
        // a) this is the root node
        // b) this is not the root and the age of the parent node is larger than the parent's age of the species node
        if ( geneNode->getAge() > min_age_left && geneNode->getAge() < max_age && geneNode->isTip() == false )
        {
            // add this node if it is within the age of our population
            nodesInPopulationSet.insert( geneNode );
        }

        if ( geneNode->isRoot() == false && ( max_age == -1.0 || max_age > geneNode->getParent().getAge() ) )
        {
            // push the parent to our current list
            individualTaxa_left.insert( &geneNode->getParent() );
        }

    }

    // get all the taxa from the species tree that are descendants of node i
    double min_age_right = n.getChild(1).getAge();
    std::vector<TopologyNode*> speciesTaxa_right;
    TreeUtilities::getTaxaInSubtree( n.getChild(1), speciesTaxa_right );

    // get all the individuals
    std::set<TopologyNode*> individualTaxa_right;
    for (size_t i = 0; i < speciesTaxa_right.size(); ++i)
    {
        const std::string &name = speciesTaxa_right[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individualTaxa_right.insert( ind[j] );
        }
    }

    // now go through all nodes in the gene
    while ( individualTaxa_right.empty() == false )
    {
        // get the first element
        std::set<TopologyNode*>::iterator it = individualTaxa_right.begin();

        // store the pointer
        TopologyNode *geneNode = *it;

        // and now remove the element from the list
        individualTaxa_right.erase( it );

        // add this node to our list of node we need to scale, if:
        // a) this is the root node
        // b) this is not the root and the age of the parent node is larger than the parent's age of the species node
        if ( geneNode->getAge() > min_age_right && geneNode->getAge() < max_age && geneNode->isTip() == false )
        {
            // add this node if it is within the age of our population
            nodesInPopulationSet.insert( geneNode );
        }

        if ( geneNode->isRoot() == false && ( max_age == -1.0 || max_age > geneNode->getParent().getAge() ) )
        {
            // push the parent to our current list
            individualTaxa_right.insert( &geneNode->getParent() );
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
void TreeNodeAgeUpdateProposal::prepareProposal( void )
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
void TreeNodeAgeUpdateProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    // no parameters

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void TreeNodeAgeUpdateProposal::undoProposal( void )
{

    // undo the proposal


    TopologyNode& parent = storedNode->getParent();

    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_new_age      = storedNode->getAge();
    double child_Age   = storedNode->getChild( 0 ).getAge();
    if ( child_Age < storedNode->getChild( 1 ).getAge())
    {
        child_Age = storedNode->getChild( 1 ).getAge();
    }

    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& geneTree = geneTrees[i]->getValue();

        std::vector<TopologyNode*> nodes = getNodesInPopulation(geneTree, *storedNode );

        for (size_t j=0; j<nodes.size(); ++j)
        {

            double new_a = nodes[j]->getAge();
            double a = new_a;
            if ( new_a > my_new_age )
            {
                a = parent_age - (parent_age - storedAge)/(parent_age - my_new_age) * (parent_age - new_a);
            }
            else
            {
                a = child_Age + (storedAge - child_Age)/(my_new_age - child_Age) * (new_a - child_Age);
            }

            // set the new age of this gene tree node
            geneTree.getNode( nodes[j]->getIndex() ).setAge( a );
        }

    }

    // set the age of the species tree node
    speciesTree->getValue().getNode( storedNode->getIndex() ).setAge( storedAge );
}


/**
 * Remove a DAG node holding a gene tree on which this move operates on.
 *
 */
void TreeNodeAgeUpdateProposal::removeGeneTree(StochasticNode<Tree> *gt)
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
void TreeNodeAgeUpdateProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


void TreeNodeAgeUpdateProposal::setProposalTuningParameter(double tp)
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
void TreeNodeAgeUpdateProposal::tune( double rate )
{

    // nothing to tune

}
