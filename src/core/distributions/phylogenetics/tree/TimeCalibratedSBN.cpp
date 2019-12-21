#include <boost/foreach.hpp>
#include "Clade.h"
#include "DistributionBeta.h"
#include "DistributionGamma.h"
#include "DistributionKumaraswamy.h"
#include "DistributionLognormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TimeCalibratedSBN.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

TimeCalibratedSBN::TimeCalibratedSBN(const SBNParameters parameters) : TypedDistribution<Tree>( new Tree() ),
    parameters( parameters ),
    taxa( parameters.getTaxa() )
{
    // Class SBNParameters handles parameterization of these edge_length_distributions
    // Here we simply use those parameters
    // Parameters are set either by calling a learn___() function or reading in an SBN
    simulateTree();

}


TimeCalibratedSBN::~TimeCalibratedSBN()
{
    // the tree will be deleted automatically by the base class

}


TimeCalibratedSBN* TimeCalibratedSBN::clone( void ) const
{

    return new TimeCalibratedSBN( *this );
}

double TimeCalibratedSBN::computeLnProbability( void )
{
    double lnProbability = 0.0;

    // Here we compute the probability of the tree topology according to the SBN

    lnProbability += parameters.computeLnProbabilityRootedTopology(*(value->clone()));

    // Add branch lengths
    lnProbability += computeLnProbabilityNodeTimes();

    return lnProbability;
}

double TimeCalibratedSBN::computeLnProbabilityNodeTimes( void )
{
    double lnProbability = 0.0;

    std::map<RbBitSet,std::vector<double> > node_age_params = parameters.getEdgeLengthDistributionParameters();

    const std::vector<TopologyNode*> tree_nodes = value->getNodes();

    // How is this SBN approximating node heights?
    bool approx_root_gamma = parameters.getBranchLengthApproximationMethod() == "rootGammaNodePropKumaraswamy";
    bool approx_prop_kumar = parameters.getBranchLengthApproximationMethod() == "rootGammaNodePropKumaraswamy";

    // We need a mapping from the ith taxon in our bitsets to the corresponding tip node age
    // Our bitsets may not match bitsets available in tree, so we make our own mapping
    std::vector<double> tip_ages;
    std::vector<Taxon> taxa = parameters.getTaxa();

    for (size_t i=0; i<taxa.size(); ++i)
    {
      size_t index = value->getTipIndex(taxa[i].getName());
      tip_ages.push_back(value->getTipNode(index).getAge());
    }

    // Get root height prob
    RbBitSet root_clade = RbBitSet(taxa.size(),true);
    std::vector<double> root_params = node_age_params[root_clade];
    double max_leaf_age = 0.0;

    for (size_t j=0; j<taxa.size(); ++j)
    {
      if ( tip_ages[j] > max_leaf_age)
      {
        max_leaf_age = tip_ages[j];
      }
    }

    if (approx_root_gamma)
    {
      double root_age_offset = value->getRoot().getAge() - max_leaf_age;
      // std::cout << "gamma shape: " << root_params.first << "; gamma rate: " << root_params.second << std::endl;
      // std::cout << "Computing probability for root age (minus offset) " << root_age_offset << ", lnProb = " << RbStatistics::Gamma::lnPdf(root_params.first, root_params.second, root_age_offset) << std::endl;
      lnProbability += RbStatistics::Gamma::lnPdf(root_params[0], root_params[1], root_age_offset);
    }
    else
    {
      throw(RbException("Invalid root age approximation scheme."));
    }

    // Get all other probs
    for (size_t i=0; i<tree_nodes.size(); ++i)
    {
      if ( !(tree_nodes[i]->isRoot() || tree_nodes[i]->isTip()) )
      {
        // Subsplit this_subsplit = tree_nodes[i]->getSubsplit(taxa);
        // RbBitSet this_split = this_subsplit.asSplitBitset();
        RbBitSet this_clade = tree_nodes[i]->getSubsplit(taxa).asCladeBitset();

        std::vector<double> these_params = node_age_params[this_clade];

        double max_descendant_leaf_age = 0.0;

        for (size_t j=0; j<taxa.size(); ++j)
        {
          if ( this_clade[j] && tip_ages[j] > max_descendant_leaf_age)
          {
            max_descendant_leaf_age = tip_ages[j];
          }
        }

        double my_age_above_tips = tree_nodes[i]->getAge() - max_descendant_leaf_age;
        double my_parents_age_above_tips = tree_nodes[i]->getParent().getAge() - max_descendant_leaf_age;
        double p = my_age_above_tips/my_parents_age_above_tips;

        if (approx_prop_kumar)
        {
          // std::cout << "Computing probability for node age with proportion p = " << p << ", lnProb = " << RbStatistics::Kumaraswamy::lnPdf(these_params[0], these_params[1], p) << std::endl;
          // std::cout << "Kumaraswamy alpha: " << these_params[0] << "; Kumaraswamy beta: " << these_params[1] << std::endl;
          lnProbability += RbStatistics::Kumaraswamy::lnPdf(these_params[0], these_params[1], p);
        }
        else
        {
          throw(RbException("Invalid node age approximation scheme."));
        }
      }
    }

    return lnProbability;
}

void TimeCalibratedSBN::redrawValue( void )
{
    simulateTree();
}

void TimeCalibratedSBN::setValue(RevBayesCore::Tree *v, bool force)
{

    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );

}


void TimeCalibratedSBN::simulateTree( void )
{

    // How is this SBN approximating node heights?
    bool approx_root_gamma = parameters.getBranchLengthApproximationMethod() == "rootGammaNodePropKumaraswamy";
    bool approx_prop_kumar = parameters.getBranchLengthApproximationMethod() == "rootGammaNodePropKumaraswamy";

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // the tree object
    // Tree *psi = value;
    Tree *psi = new Tree();

    psi->setRooted(true);

    // For drawing branch lengths
    std::map<RbBitSet,std::vector<double> > node_age_params = parameters.getNodeTimeDistributionParameters();

    // create the tip nodes, this will always result in tip node i being the ith taxon in our ordered set (convenient)
    std::vector<TopologyNode*> tip_nodes;
    std::vector<double> tip_ages;
    double max_leaf_age = 0.0;
    for (size_t i=0; i<taxa.size(); ++i)
    {
        // Node age
        double node_age = taxa[i].getAge();
        tip_ages.push_back(node_age);

        if (node_age > max_leaf_age)
        {
          max_leaf_age = node_age;
        }

        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );
        node->setAge(node_age);
        tip_nodes.push_back(node);
    }

    // List of active tree nodes/subsplits
    // We pair them such that each tree node corresponds to the subsplit it defines
    std::vector<std::pair<Subsplit,TopologyNode*> > active;

    // Root split
    double u = rng->uniform01();
    TopologyNode* root = new TopologyNode();
    root->setNodeType(false, true, false);
    Subsplit root_split = parameters.drawRootSplit();
    active.push_back(std::make_pair(root_split,root));

    // Draw root age
    RbBitSet root_clade = RbBitSet(taxa.size(),true);
    std::vector<double> root_params = node_age_params[root_clade];

    if ( approx_root_gamma )
    {
      double offset_age  = RbStatistics::Gamma::rv(root_params[0], root_params[1], *rng);
      root->setAge(offset_age + max_leaf_age);
    }
    else
    {
      throw(RbException("Invalid root height approximation."));
    }

    // All other subplits
    while (active.size() > 0)
    {
      // Get a node/subsplit to work on, remove that from list
      std::pair<Subsplit,TopologyNode*> this_parent = active.back();
      active.pop_back();

      Subsplit this_parent_subsplit = this_parent.first;
      TopologyNode* this_parent_node = this_parent.second;

      TopologyNode* Y_child_node;
      TopologyNode* Z_child_node;

      // Choose subsplit of Y
      Subsplit Y_child = parameters.drawSubsplitForY(this_parent_subsplit);
      if ( Y_child.isFake() )
      {
        // This is a tip, we don't add it to the active pile, nor do we give it an age
        Y_child_node = tip_nodes[Y_child.getFsbY()];
      }
      else
      {
        // This is an internal node
        Y_child_node = new TopologyNode();
        active.push_back(std::make_pair(Y_child,Y_child_node));

        RbBitSet y_clade = Y_child.asCladeBitset();

        std::vector<double> these_params = node_age_params[y_clade];

        double max_descendant_leaf_age = 0.0;

        for (size_t j=0; j<taxa.size(); ++j)
        {
          if ( y_clade[j] && tip_ages[j] > max_descendant_leaf_age)
          {
            max_descendant_leaf_age = tip_ages[j];
          }
        }

        double p;
        if ( approx_prop_kumar )
        {
          p  = RbStatistics::Kumaraswamy::rv(these_params[0], these_params[1], *rng);
        }
        else
        {
          throw(RbException("Invalid node age approximation scheme."));
        }

        double my_parents_age_above_tips = this_parent_node->getAge();

        double my_age_above_tips = max_descendant_leaf_age + p * (my_parents_age_above_tips - max_descendant_leaf_age);

        Y_child_node->setAge(my_age_above_tips);

      }

      // Choose subsplit of Z
      Subsplit Z_child = parameters.drawSubsplitForZ(this_parent_subsplit);
      if ( Z_child.isFake() )
      {
        // This is a tip, we don't add it to the active pile, nor do we give it an age
        Z_child_node = tip_nodes[Z_child.getFsbY()];
      }
      else
      {
        // This is an internal node
        Z_child_node = new TopologyNode();
        active.push_back(std::make_pair(Z_child,Z_child_node));

        RbBitSet z_clade = Z_child.asCladeBitset();

        std::vector<double> these_params = node_age_params[z_clade];

        double max_descendant_leaf_age = 0.0;

        for (size_t j=0; j<taxa.size(); ++j)
        {
          if ( z_clade[j] && tip_ages[j] > max_descendant_leaf_age)
          {
            max_descendant_leaf_age = tip_ages[j];
          }
        }

        double p;
        if ( approx_prop_kumar )
        {
          p  = RbStatistics::Kumaraswamy::rv(these_params[0], these_params[1], *rng);
        }
        else
        {
          throw(RbException("Invalid node age approximation scheme."));
        }

        double my_parents_age_above_tips = this_parent_node->getAge();

        double my_age_above_tips = max_descendant_leaf_age + p * (my_parents_age_above_tips - max_descendant_leaf_age);

        Z_child_node->setAge(my_age_above_tips);

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

    psi->orderNodesByIndex();

  // finally store the new value
    delete value;
    value = psi;

}


/** Swap a parameter of the distribution */
void TimeCalibratedSBN::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
  // SBNs have their parameters set once and they do not change afterwards
}

// /**
//  * We check here if all the constraints are satisfied.
//  * These are hard constraints, that is, the clades must be monophyletic.
//  *
//  * \return     True if the constraints are matched, false otherwise.
//  */
// bool TimeCalibratedSBN::matchesConstraints( void )
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
