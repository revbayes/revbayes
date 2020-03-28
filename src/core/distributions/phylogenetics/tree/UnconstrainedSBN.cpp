#include <boost/foreach.hpp>
#include "Clade.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "TopologyNode.h"
#include "UnconstrainedSBN.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

// Empty constructor (needed by smart moves)
UnconstrainedSBN::UnconstrainedSBN(void) : TypedDistribution<Tree>( new Tree() ),
    parameters( ),
    taxa( )
{

}

UnconstrainedSBN::UnconstrainedSBN(const SBNParameters parameters) : TypedDistribution<Tree>( new Tree() ),
    parameters( parameters ),
    taxa( parameters.getTaxa() )
{
    // Class SBNParameters handles parameterization of these edge_length_distributions
    // Here we simply use those parameters
    // Parameters are set either by calling a learn___() function or reading in an SBN
    simulateTree();

}


UnconstrainedSBN::~UnconstrainedSBN()
{
    // the tree will be deleted automatically by the base class

}


UnconstrainedSBN* UnconstrainedSBN::clone( void ) const
{

    return new UnconstrainedSBN( *this );
}

double UnconstrainedSBN::computeLnProbability( void )
{
    double lnProbability = 0.0;

    // Here we compute the probability of the tree topology according to the SBN

    // lnProbability += computeLnProbabilityUnrootedTopologyMarginalize();
    lnProbability += parameters.computeLnProbabilityUnrootedTopology(*value);

    // Add branch lengths
    lnProbability += computeLnProbabilityBranchLengths();

    return lnProbability;
}

double UnconstrainedSBN::computeLnProbabilityBranchLengths( void )
{
    double lnProbability = 0.0;

    // Get branch lengths
    const std::vector<TopologyNode*> tree_nodes = value->getNodes();
    for (size_t i=0; i<tree_nodes.size(); ++i)
    {
      if (!tree_nodes[i]->isRoot())
      {
        RbBitSet this_split = tree_nodes[i]->getSubsplit(taxa).asSplitBitset();
        lnProbability += parameters.computeEdgeLengthProbability(this_split,tree_nodes[i]->getBranchLength());
      }
    }

    return lnProbability;
}

double UnconstrainedSBN::computeLnProbabilityUnrootedTopologyMarginalize( void )
{
    // Make the tree properly rooted so that rooting to branches has the desired effect

    double lnProbability = 0.0;

    std::vector<double> lnl_given_root;
    double offset = RbConstants::Double::neginf;

    // sum over rooting locations
    for (size_t ri=0; ri < value->getNumberOfNodes(); ++ri)
    {
      if (!(value->getNode(ri).isRoot()))
      {
        Tree *tree_rooted = value->clone();
        tree_rooted->makeRooted(tree_rooted->getNode(ri),true);
        double lnl = parameters.computeLnProbabilityRootedTopology( *(tree_rooted) );
        lnl_given_root.push_back(lnl);
        delete tree_rooted;
      }
    }

    lnProbability = RbMath::log_sum_exp(lnl_given_root);

    return lnProbability;
}

void UnconstrainedSBN::redrawValue( void )
{
    simulateTree();
}

void UnconstrainedSBN::setValue(RevBayesCore::Tree *v, bool force)
{

    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );

}


void UnconstrainedSBN::simulateTree( void )
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // the tree object
    // Tree *psi = value;
    Tree *psi = new Tree();

    // We always draw a rooted tree, if we want to unroot we need it to first be rooted
    psi->setRooted(true);

    // create the tip nodes
    std::vector<TopologyNode*> tip_nodes;
    for (size_t i=0; i<taxa.size(); ++i)
    {
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );
        tip_nodes.push_back(node);
    }

    // List of active tree nodes/subsplits
    // We pair them such that each tree node corresponds to the subsplit it defines
    std::vector<std::pair<size_t,TopologyNode*> > active;

    // Root split
    double u = rng->uniform01();
    TopologyNode* root = new TopologyNode();
    root->setNodeType(false, true, false);
    size_t root_split = parameters.drawRootSplit();
    active.push_back(std::make_pair(root_split,root));

    // All other subplits
    // TODO: we should exploit the fact that tips are in the master list and have 0 children to avoid "getting" subsplits just to check if they're tips
    while (active.size() > 0)
    {
      // Get a node/subsplit to work on, remove that from list
      std::pair<size_t,TopologyNode*> this_parent = active.back();
      active.pop_back();

      size_t this_parent_subsplit = this_parent.first;
      TopologyNode* this_parent_node = this_parent.second;

      TopologyNode* Y_child_node;
      TopologyNode* Z_child_node;

      // Choose subsplit of Y
      size_t Y_child = parameters.drawSubsplitForY(this_parent_subsplit);
      if ( parameters.getNumChildrenForParent(Y_child) == 0 )
      {
        // This is a tip, we don't add it to the active pile
        Y_child_node = tip_nodes[parameters.getSubsplitReference(Y_child).getFsbY()];
      }
      else
      {
        // This is an internal node
        Y_child_node = new TopologyNode();
        active.push_back(std::make_pair(Y_child,Y_child_node));
      }

      // Choose subsplit of Z
      size_t Z_child = parameters.drawSubsplitForZ(this_parent_subsplit);
      if ( parameters.getNumChildrenForParent(Z_child) == 0 )
      {
        // This is a tip, we don't add it to the active pile
        Z_child_node = tip_nodes[parameters.getSubsplitReference(Z_child).getFsbY()];
      }
      else
      {
        // This is an internal node
        Z_child_node = new TopologyNode();
        active.push_back(std::make_pair(Z_child,Z_child_node));
      }

      // Attach nodes to eachother
      this_parent_node->addChild(Y_child_node);
      this_parent_node->addChild(Z_child_node);
      Y_child_node->setParent(this_parent_node);
      Z_child_node->setParent(this_parent_node);
    }

    // initialize the topology by setting the root
    psi->setRoot(root, true);

    // TODO: with the selection of taxa for tips using bitsets, taxon names are paired to indices, the following may be extraneous
    // re-couple tip node names with tip indices
    // this is necessary because otherwise tip names get scrambled across replicates
    for (size_t i=0; i<taxa.size(); i++)
    {
    	psi->getTipNodeWithName(taxa[i].getName()).setIndex(i);
    }

    psi->unroot();
    psi->setRooted( false );
    psi->orderNodesByIndex();

    // Add branch lengths
    const std::vector<TopologyNode*> tree_nodes = psi->getNodes();
    for (size_t i=0; i<tree_nodes.size(); ++i)
    {
      if (!tree_nodes[i]->isRoot())
      {
        // Subsplit this_subsplit = tree_nodes[i]->getSubsplit(taxa);
        // RbBitSet this_split = this_subsplit.asSplitBitset();
        RbBitSet this_split = tree_nodes[i]->getSubsplit(taxa).asSplitBitset();
        double brlen = parameters.drawEdgeLength(this_split);
        tree_nodes[i]->setBranchLength(brlen,false);
      }
    }

  // finally store the new value
    delete value;
    value = psi;

}


/** Swap a parameter of the distribution */
void UnconstrainedSBN::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
  // SBNs have their parameters set once and they do not change afterwards
}

// /**
//  * We check here if all the constraints are satisfied.
//  * These are hard constraints, that is, the clades must be monophyletic.
//  *
//  * \return     True if the constraints are matched, false otherwise.
//  */
// bool UnconstrainedSBN::matchesConstraints( void )
// {
//
//     if ( constraints.empty() == true )
//     {
// 		return true;
// 	}
//     else
//     {
//
// 		const TopologyNode &root = value->getRoot();
// 		for (std::vector<Clade>::iterator it = constraints.begin(); it != constraints.end(); ++it)
// 		{
// 			if ( root.containsClade( *it, true ) == false )
// 			{
// 				return false;
// 			}
// 		}
//
// 		return true;
// 	}
//
// }
