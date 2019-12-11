/**
 * @file
 * This file contains the declaration of SBNParameters, which is
 * class that holds the parameters for an SBN, and methods for use
 * with SBNs.
 *
 * @brief Implementation of SBNParameters
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#include <boost/foreach.hpp>
#include "DistributionGamma.h"
#include "DistributionGeneralizedGamma.h"
#include "DistributionKumaraswamy.h"
#include "DistributionLognormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "Split.h"
#include "SBNParameters.h"

#include <cmath>
#include <iomanip>


using namespace RevBayesCore;


/** Construct empty SBN parameters */
SBNParameters::SBNParameters( void ) :
  branch_length_approximation_method(),
  edge_length_distribution_parameters(),
  num_taxa(),
  root_splits(),
  subsplit_cpds(),
  taxa(),
  time_calibrated()
{
  // Nothing to do
}

/** Construct empty SBN parameters from taxa */
SBNParameters::SBNParameters( std::vector<Taxon> taxa, const std::string &branch_length_approximation, bool allow_unseen, double epsilon ) :
  branch_length_approximation_method( branch_length_approximation ),
  edge_length_distribution_parameters(),
  num_taxa( taxa.size() ),
  root_splits(),
  subsplit_cpds(),
  taxa( taxa ),
  time_calibrated()
{

}

/** Construct rate matrix with n states */
SBNParameters::SBNParameters( const SBNParameters &sbn ) :
  branch_length_approximation_method(sbn.branch_length_approximation_method),
  edge_length_distribution_parameters(sbn.edge_length_distribution_parameters),
  num_taxa(sbn.num_taxa),
  time_calibrated(sbn.time_calibrated),
  root_splits(sbn.root_splits),
  taxa(sbn.taxa),
  subsplit_cpds(sbn.subsplit_cpds)
{

}

SBNParameters::~SBNParameters()
{

}


/** Construct rate matrix with n states */
SBNParameters& SBNParameters::operator=( const SBNParameters &sbn )
{

    if ( this != &sbn )
    {
      branch_length_approximation_method  = sbn.branch_length_approximation_method;
      edge_length_distribution_parameters = sbn.edge_length_distribution_parameters;
      num_taxa                            = sbn.num_taxa;
      root_splits                         = sbn.root_splits;
      subsplit_cpds                       = sbn.subsplit_cpds;
      taxa                                = sbn.taxa;
      time_calibrated                     = sbn.time_calibrated;
    }

    return *this;
}

const std::string& SBNParameters::getBranchLengthApproximationMethod(void) const
{
  return branch_length_approximation_method;
}

std::map<RbBitSet,std::vector<double> >& SBNParameters::getEdgeLengthDistributionParameters(void)
{
  return edge_length_distribution_parameters;
}

std::map<RbBitSet,std::vector<double> >& SBNParameters::getNodeTimeDistributionParameters(void)
{
  return edge_length_distribution_parameters;
}

const size_t SBNParameters::getNumTaxa(void) const
{
  return num_taxa;
}

std::vector<std::pair<Subsplit,double> >& SBNParameters::getRootSplits(void)
{
  return root_splits;
}

std::map<Subsplit,std::vector<std::pair<Subsplit,double> > >& SBNParameters::getSubsplitCPDs(void)
{
  return subsplit_cpds;
}

std::vector<Taxon>& SBNParameters::getTaxa(void)
{
  return taxa;
}

const std::vector<Taxon>& SBNParameters::getTaxa(void) const
{
  return taxa;
}

// Needed for fitting a number of branch-length distributions
std::vector<double> get_moments(std::vector<double> samples)
{
  // Get mean/var of observations
  double mean;
  for (size_t i=0; i<samples.size(); ++i)
  {
    mean += samples[i];
  }
  mean /= samples.size();

  double var;
  for (size_t i=0; i<samples.size(); ++i)
  {
    var += pow(samples[i] - mean,2.0);
  }
  var /= samples.size();

  std::vector<double> moments;
  moments.push_back(mean);
  moments.push_back(var);

  return moments;
}

//Need for fitting time-SBNs and unrooted SBNs
std::vector<double> fit_gamma_MOM(std::vector<double> samples)
{
  std::vector<double> moments = std::vector<double>(2,-1.0);
  moments = get_moments(samples);
  std::vector<double> par = std::vector<double>(2,-1.0);
  par[1] = moments[0]/moments[1];
  par[0] = moments[0]*moments[0]/moments[1];


// std::cout << "fitting gamma via MOM" << std::endl;
// std::cout << "moments are mean = " << moments[0] << " and var = " << moments[1] << std::endl;
// std::cout << "params are alpha = " << par[0] << " and beta = " << par[1] << std::endl;

  return par;
}

// Need for fitting time-SBNs
// Jones (2009) eqn 5.2, Score function for Kumaraswamy distribution
double kumar_score(double a, std::vector<double> &samples)
{
  double n = (double)samples.size();

  double T1 = 0.0;
  double T2 = 0.0;
  double T3 = 0.0;

  for (size_t i=0; i<samples.size(); ++i)
  {
    double y = pow(samples[i],a);
    double log_y = log(y);
    double one_minus_y = 1.0 - y;
    T1 += log_y/one_minus_y;
    T2 += (y * log_y)/one_minus_y;
    T3 += log(one_minus_y);
  }

  T1 *= 1.0/n;
  T2 *= 1.0/n;
  T3 *= 1.0/n;

  return n/a * (1 + T1 + T2/T3);

}

// Need for fitting time-SBNs
// Fits a Kumaraswamy distribution via Nesterov's accelerated gradient descent
std::vector<double> fit_kumar_agd(std::vector<double> &samples)
{

  double x_s = 1.0;
  double x_s_1;

  double y_s = 1.0;
  double y_s_1;

  double lambda_s = 1.0;
  double lambda_s_1 = (1.0 + sqrt(1.0 + 4.0 * lambda_s*lambda_s)) / 2.0;

  double d = 10000000.0;
  double beta_inv = 1.0/1000.0;
  double eps = 0.001;
  // Run accelerated gradient descent
  while (d > eps)
  {
    double gamma_s = (1.0 - lambda_s) / lambda_s_1;

    y_s_1 = x_s + beta_inv * kumar_score(x_s,samples);
    x_s_1 = (1.0 - gamma_s) * y_s_1 + gamma_s * y_s;

    d = fabs(x_s_1 - x_s);

    x_s = x_s_1;
    y_s = y_s_1;

    lambda_s = lambda_s_1;
    lambda_s_1 = (1.0 + sqrt(1.0 + 4.0 * lambda_s*lambda_s)) / 2.0;

  }

  std::vector<double> alphabeta;
  alphabeta.push_back(x_s_1);

  // Get beta parameter of Kumaraswamy distribution (x_s_1 is alpha)
  double s = 0.0;
  for (size_t i=0; i<samples.size(); ++i)
  {
    s += log(1.0 - pow(samples[i],x_s_1));
  }

  alphabeta.push_back(-(double)samples.size()/s);

  return alphabeta;
}

// Need for unrooted SBNs to fit generalized gamma distributions
double fn_fixed_c_gengamma_lnl(double c, std::vector<double> samples)
{
// std::cout << "calling fn_fixed_c_gengamma_lnl with c = " << c << std::endl;
  // y = x^c
  std::vector<double> y;
  for (size_t i=0; i<samples.size(); ++i)
  {
    y.push_back(pow(samples[i],c));
  }
  // ypar = gammaMoments2Params(mean(y),var(y))
  std::vector<double> transformed_par = std::vector<double>(2,-1.0);
  transformed_par = fit_gamma_MOM(y);
// std::cout << "transformed_par[0] = " << transformed_par[0] << "; transformed_par[1] = " << transformed_par[1] << std::endl;
  double l = transformed_par[0];
  double a = pow(transformed_par[1],(1.0/c));

// std::cout << "a = " << a << "; l = " << l << std::endl;
  double lnL = 0.0;
  for (size_t i=0; i<samples.size(); ++i)
  {
    lnL += RbStatistics::GeneralizedGamma::lnPdf(a,c,l,samples[i]);
  }
// std::cout << "fn_fixed_c_gengamma_lnl is " << lnL << " for a = " << a << "; c = " << c << "; l = " << l << std::endl;
  return lnL;
}

// Need for unrooted SBNs
// Fits a Generalized Gamma distribution by optimizing over the c parameter, for any given value fits a,l by MOM
std::vector<double> fit_gengam_gss(std::vector<double> samples)
{

  // Constants
  double eps = 1e-5; // tolerance for minimum
  double gold_ratio = (1.0 + sqrt(5.0))/2.0; // the golden ratio (should eventually be moved to RbConstants)

  double x1 = 0.01;
  double x2 = 5.0;

  double f_x1 = -fn_fixed_c_gengamma_lnl(x1,samples);
  double f_x2 = -fn_fixed_c_gengamma_lnl(x2,samples);

  double x3 = x2 - (x2 - x1)/gold_ratio;
  double x4 = x1 + (x2 - x1)/gold_ratio;

  double f_x3 = -fn_fixed_c_gengamma_lnl(x3,samples);
  double f_x4 = -fn_fixed_c_gengamma_lnl(x4,samples);

// std::cout << "starting search, f_x1 = " << f_x1 << "; f_x2 = " << f_x2 << "; f_x3 = " << f_x3  << "; f_x4 = " << f_x4 << std::endl;
  while ( fabs(x3 - x4) > eps)
  {
    if (f_x3 < f_x4) {
      x2 = x4;
      x4 = x3;
      f_x4 = f_x3;
      x3 = x2 - (x2 - x1)/gold_ratio;
      x4 = x1 + (x2 - x1)/gold_ratio;
      f_x3 = -fn_fixed_c_gengamma_lnl(x3,samples);
    } else {
      x1 = x3;
      x3 = x4;
      f_x3 = f_x4;
      x3 = x2 - (x2 - x1)/gold_ratio;
      x4 = x1 + (x2 - x1)/gold_ratio;
      f_x4 = -fn_fixed_c_gengamma_lnl(x4,samples);
    }
// std::cout << "    x3 = " << x3 << "; x4 = " << x4 << "; f_x3 = " << f_x3 << "; f_x4 = " << f_x4 << std::endl;
  }

  // Best value
  double c = (x1+x2)/2.0;

  // Get a,l
  std::vector<double> y;
  for (size_t i=0; i<samples.size(); ++i)
  {
    y.push_back(pow(samples[i],c));
  }

  std::vector<double> transformed_par = std::vector<double>(2,-1.0);
  transformed_par = fit_gamma_MOM(y);
  double l = transformed_par[0];
  double a = pow(transformed_par[1],(1.0/c));

  std::vector<double> par;
  par.push_back(a);
  par.push_back(c);
  par.push_back(l);
  // std::cout << "completed fit_gengam_gss, got parameters par[0] = " << par[0] << "; par[1] = " << par[1]  << "; par[2] = " << par[2] << std::endl;
  return par;
}


/* Computes the probability of the given ROOTED tree under this SBN. */
double SBNParameters::computeLnProbabilityRootedTopology( const Tree &tree ) const
{
  if ( !(tree.isBinary()) )
  {
    std::cout << "tree no longer binary" << std::endl;
    return RbConstants::Double::nan;
  }

  double lnProbability = 0.0;

  // first compute root split probability
  const TopologyNode &root = tree.getRoot();
  const std::vector<TopologyNode*>& root_children = root.getChildren();

  if ( root_children.size() != 2 )
  {
    return( RbConstants::Double::nan );
  }

  Subsplit root_split = tree.getRootSubsplit(taxa);

  lnProbability += computeRootSplitProbability(root_split);

  // get all parent-child subsplit pairs, calculate their probabilities
  // TODO: We could do this more efficiently by travesing the tree, since we would not have to find the subsplit every node belongs to every time
  std::vector<std::pair<Subsplit,Subsplit> > parent_child_subsplits = tree.getAllSubsplitParentChildPairs(taxa);

  std::pair<Subsplit,Subsplit> parent_child_pair;

  BOOST_FOREACH(parent_child_pair, parent_child_subsplits)
  {
    lnProbability += computeSubsplitTransitionProbability(parent_child_pair.first, parent_child_pair.second);
  }

  return lnProbability;

}

/* Helper function for UNROOTED trees. Computes for each branch the joint probability Pr(tree and root here) */
std::vector<std::pair<Subsplit,double> > SBNParameters::computeLnProbabilityTopologyAndRooting( const Tree &t ) const
{
    std::vector<std::pair<Subsplit,double> > lnProbability;

    Tree* tree = t.clone();

    // Prep for tip to root pass
    std::string order = "postorder";
    tree->orderNodesForTraversal(order);
    const std::vector<TopologyNode*> &postorder_nodes = tree->getNodes();

    // For storing subsplits
    std::vector<Subsplit> per_node_subsplit = std::vector<Subsplit>(tree->getNumberOfNodes(),Subsplit());
    // For storing cumulative log-probabilities from tip through this branch.
    std::vector<double> ttr = std::vector<double>(tree->getNumberOfNodes(),0.0);

    // Tip to root pass, here we do two things
    // 1) Get all nodes' subsplits (we will need these repeatedly)
    //      We do not make a subsplit for the root as the root is a trifurcation in an unrooted tree
    // 2) accumulate tipward tree probabilities, ttr[node] = log(Pr(all subsplits tipward of this subsplit | this subsplit))
    for (std::vector<TopologyNode*>::const_iterator it = postorder_nodes.begin(); it != (postorder_nodes.end()-1); ++it)
    {
      size_t index = (*it)->getIndex();

      if ( (*it)->isTip() )
      {
        // 1)
        per_node_subsplit[index] = (*it)->getSubsplit(taxa);

        // 2)
        ttr[index] = 0.0;
      }
      else
      {
        std::vector<int> children = (*it)->getChildrenIndices();

        // 1)
        RbBitSet clade_1 = per_node_subsplit[children[0]].asCladeBitset();
        RbBitSet clade_2 = per_node_subsplit[children[1]].asCladeBitset();
        per_node_subsplit[index] = Subsplit(clade_1,clade_2);

        // 2)
        ttr[index] = computeSubsplitTransitionProbability(per_node_subsplit[index],per_node_subsplit[children[0]]) + ttr[children[0]] + computeSubsplitTransitionProbability(per_node_subsplit[index],per_node_subsplit[children[1]]) + ttr[children[1]];
      }
    }

    // Root to tip pass (this is where the fun starts)
    // Here we do several things
    // 1) We collect a vector at every node, rtt[node], which handles cumulative tree probabilities like ttr[node], but for the complement of the tree tipward of node.
    //    rtt[node] is a log-scale probability
    //    If node is a root descendant, we start the propagation.
    // 2) We use rtt[node], ttr[node], and Pr(root on this edge) to compute Pr(tree and root here)
    //
    order = "preorder";
    tree->orderNodesForTraversal(order);
    const std::vector<TopologyNode*> &preorder_nodes = tree->getNodes();

    std::vector<double> rtt = std::vector<double>(tree->getNumberOfNodes(),0.0);

    // // For storing subsplits from alternative orientations of the tree
    // std::vector<Subsplit> per_node_reversed_subsplit = std::vector<Subsplit>(tree->getNumberOfNodes(),Subsplit());

    // Loop over edges of tree (exploit equivalency between an edge and the node that edge subtends)
    // The root has no edge so there is nothing to do for the root, so we skip it
    for (std::vector<TopologyNode*>::const_iterator it = preorder_nodes.begin()+1; it != preorder_nodes.end(); ++it)
    {
      // std::cout << ">>>working on a root/internal/tip node " << ((*it)->isRoot()) << "/" << ((*it)->isInternal()) << "/" << ((*it)->isTip()) << std::endl;
      // std::cout << ">The node's subsplit is " << per_node_subsplit[(*it)->getIndex()] << std::endl;
      // std::cout << ">The node's parent's subsplit is " << per_node_subsplit[(*it)->getParent().getIndex()] << std::endl;

      size_t index = (*it)->getIndex();

      // Edges descending from root need to be handled differently
      if ( (*it)->getParent().isRoot() )
      {
        // Get subsplits for other two descendants of root
        std::vector<int> root_children_indices = tree->getRoot().getChildrenIndices();

        std::vector<int> sibling_indices;

        for (size_t i=0; i<3; ++i)
        {
          if (index != root_children_indices[i])
          {
            sibling_indices.push_back(root_children_indices[i]);
          }
        }

        RbBitSet sib_0_bitset = per_node_subsplit[sibling_indices[0]].asCladeBitset();
        RbBitSet sib_1_bitset = per_node_subsplit[sibling_indices[1]].asCladeBitset();
        Subsplit sib_split = Subsplit(sib_0_bitset,sib_1_bitset);

        std::vector<std::pair<Subsplit,Subsplit> > cases = per_node_subsplit[index].doVirtualRootingRootParent(per_node_subsplit[sibling_indices[0]],per_node_subsplit[sibling_indices[1]],per_node_subsplit[(*it)->getIndex()]);

        // 1) Propagate rtt
        rtt[index] = ttr[sibling_indices[0]] + computeSubsplitTransitionProbability(sib_split,per_node_subsplit[sibling_indices[0]]) + ttr[sibling_indices[1]] + computeSubsplitTransitionProbability(sib_split,per_node_subsplit[sibling_indices[1]]);

        // 2) Compute Pr(tree,root here)
        double pr_subtree_s = computeSubsplitTransitionProbability(cases[2].first,cases[2].second) + ttr[index];
        double pr_subtree_not_s = computeSubsplitTransitionProbability(cases[5].first,cases[5].second) + rtt[index];
        std::pair<Subsplit,double> root_split_and_tree_prob;
        root_split_and_tree_prob.first = cases[2].first;
        root_split_and_tree_prob.second = pr_subtree_s + pr_subtree_not_s + computeRootSplitProbability(cases[2].first);
        lnProbability.push_back(root_split_and_tree_prob);

      }
      else
      {
        // std::cout << "(*it)->getParent().isRoot() == " << (*it)->getParent().isRoot() << std::endl;
        // Define parent-child pair for current rooting (parent first, child second)
        // std::pair<Subsplit,Subsplit> this_parent_child;
        // this_parent_child.first = per_node_subsplit[(*it)->getParent().getIndex()];
        // this_parent_child.second = per_node_subsplit[(*it)->getParent().getIndex()];

        // Get all cases for virtual rooting of this edge (including current rooting)
        // std::vector<std::pair<Subsplit,Subsplit> > cases = per_node_subsplit[index].doVirtualRootingNonRootParent(this_parent_child.first,this_parent_child.second);
        std::vector<std::pair<Subsplit,Subsplit> > cases = per_node_subsplit[index].doVirtualRootingNonRootParent(per_node_subsplit[(*it)->getParent().getIndex()],per_node_subsplit[(*it)->getIndex()]);

        // for (size_t i=0; i<6; ++i) {
        //   std::cout << "cases[" << i << "].first = " << cases[i].first << std::endl;
        //   std::cout << "cases[" << i << "].second = " << cases[i].second << std::endl;
        // }

        // We will need to get the first subsplit both in our parent's sister node and in the rest of the tree
        // Thus we need our sibling's index and our parent's sibling's index
        size_t parent_index = (*it)->getParent().getIndex();
        size_t grandparent_index = (*it)->getParent().getParent().getIndex();
        size_t my_siblings_index = 0;
        size_t my_parents_siblings_index = 0;

        if ( &((*it)->getParent().getChild(0)) == (*it) )
        {
          my_siblings_index = 1;
        }
        my_siblings_index = (*it)->getParent().getChild(my_siblings_index).getIndex();

        if ( &((*it)->getParent().getParent().getChild(0)) == &((*it)->getParent()) )
        {
          my_parents_siblings_index = 1;
        }
        my_parents_siblings_index = (*it)->getParent().getParent().getChild(my_parents_siblings_index).getIndex();

        // 1) Propagate rtt
        RbBitSet parents_clade = per_node_subsplit[(*it)->getParent().getIndex()].asCladeBitset();
        RbBitSet parents_sisters_clade = per_node_subsplit[my_parents_siblings_index].asCladeBitset();
        RbBitSet not_parents_clade = parents_clade;
        ~not_parents_clade;
        RbBitSet not_parents_sisters_clade = parents_sisters_clade;
        ~not_parents_sisters_clade;
        RbBitSet everyone_else = not_parents_clade & not_parents_sisters_clade;
        Subsplit everyone_else_from_parents_sister = Subsplit(everyone_else,parents_sisters_clade);
// std::cout << "me             = " << per_node_subsplit[index].asCladeBitset() << std::endl;
// std::cout << "parent         = " << parents_clade << std::endl;
// std::cout << "parents_sister = " << parents_sisters_clade << std::endl;
// std::cout << "everyone_else  = " << everyone_else << std::endl;
// std::cout << "ee || ps       = " << everyone_else_from_parents_sister << std::endl;
//
// std::cout << "rtt[parent_index] = " << rtt[parent_index] << std::endl;
// std::cout << "looking for subsplit transition probability for parent->child split " << cases[5].second << " -> " << everyone_else_from_parents_sister << std::endl;
// std::cout << "computeSubsplitTransitionProbability(cases[5].second,everyone_else_from_parents_sister) = " << computeSubsplitTransitionProbability(cases[5].second,everyone_else_from_parents_sister) << std::endl;
// std::cout << "ttr[my_siblings_index] = " << ttr[my_siblings_index] << std::endl;
// std::cout << "looking for subsplit transition probability for parent->child split " << cases[5].second << " -> " << per_node_subsplit[my_siblings_index] << std::endl;
// std::cout << "computeSubsplitTransitionProbability(cases[5].second,per_node_subsplit[my_siblings_index]) = " << computeSubsplitTransitionProbability(cases[5].second,per_node_subsplit[my_siblings_index]) << std::endl;
        rtt[index] = ttr[my_siblings_index] + computeSubsplitTransitionProbability(cases[5].second,everyone_else_from_parents_sister) + rtt[parent_index] + computeSubsplitTransitionProbability(cases[5].second,per_node_subsplit[my_siblings_index]);

        // 2) Compute Pr(tree,root here)
        double pr_subtree_s = computeSubsplitTransitionProbability(cases[2].first,cases[2].second) + ttr[index];
        double pr_subtree_not_s = computeSubsplitTransitionProbability(cases[5].first,cases[5].second) + rtt[index];
        std::pair<Subsplit,double> root_split_and_tree_prob;
        root_split_and_tree_prob.first = cases[2].first;
        root_split_and_tree_prob.second = pr_subtree_s + pr_subtree_not_s + computeRootSplitProbability(cases[2].first);
        lnProbability.push_back(root_split_and_tree_prob);
      }

    }
    delete tree;

  return lnProbability;
}

/* Computes the probability of the given UNROOTED tree under this SBN. */
double SBNParameters::computeLnProbabilityUnrootedTopology( const Tree &tree ) const
{
  std::vector<std::pair<Subsplit,double> > lnPr_tree_and_rooting = computeLnProbabilityTopologyAndRooting(tree);

  double lnProbability = 0.0;

  std::vector<double> per_edge_log_probs;
  std::pair<Subsplit,double> this_rooting;
  // Loop over parent subsplits
std::cout << "per_edge_log_probs = ";
  BOOST_FOREACH(this_rooting, lnPr_tree_and_rooting) {
    per_edge_log_probs.push_back(this_rooting.second);
std::cout << this_rooting.second << ", ";
  }
std::cout << std::endl;

  lnProbability = RbMath::log_sum_exp(per_edge_log_probs);

  return lnProbability;
}

/* Computes the probability of seeing a particular root split given an SBN */
double SBNParameters::computeRootSplitProbability( const Subsplit &root_split ) const
{
  double log_prob = RbConstants::Double::neginf;
  for (size_t i=0; i<root_splits.size(); ++i)
  {
    if ( root_split == root_splits[i].first)
    {
      log_prob = log(root_splits[i].second);
    }
  }
  return log_prob;
}

/* Computes the probability of seeing a particular parent-child subsplit pair given an SBN */
double SBNParameters::computeSubsplitTransitionProbability( const Subsplit &parent, const Subsplit &child ) const
{

  // Find all potential children of parent
  const std::vector<std::pair<Subsplit,double> > &all_children = subsplit_cpds.at(parent);

  double log_prob = RbConstants::Double::neginf;

  for (size_t i=0; i<all_children.size(); ++i)
  {
    if ( child == all_children[i].first)
    {
      log_prob = log(all_children[i].second);
    }
  }
  return log_prob;

}

// /* Computes the probability of seeing a particular parent-child subsplit pair given an SBN */
// std::vector<std::pair<Split,double> > SBNParameters::computeSplitProbabilities( void ) const
// {
//
//   // Find all potential children of parent
//   std::vector<std::pair<Split,double> > split_probs;
//
//   return split_probs;
//
// }

/* Draw a root split from an SBN */
Subsplit SBNParameters::drawRootSplit( void ) const
{

  double u = GLOBAL_RNG->uniform01();
  size_t index;
  for (size_t i=0; i<root_splits.size(); ++i)
  {
    if (u < root_splits[i].second)
    {
      index = i;
      break;
    }
    u -= root_splits[i].second;
  }

  return root_splits[index].first;
}

/* Draws a subsplit for subsplit S's clade Y*/
Subsplit SBNParameters::drawSubsplitForY( const Subsplit &s ) const
{

  // Find all potential children of s
  std::vector<std::pair<Subsplit,double> > my_children = subsplit_cpds.at(s);

  // Find a distinguishing feature of clade Y in subsplit s
  // Since Y and Z are disjoint, we can use the first set bit in Y
  size_t fsb = s.getYBitset().getFirstSetBit();

  double u = GLOBAL_RNG->uniform01();
  size_t index;

  // This is like drawing a root split, but we must ensure we only draw from subsplits of Y
  for (size_t i=0; i<my_children.size(); ++i)
  {
    // This is a subsplit of Y if one of its splits has the same first set bit as Y
    // my_children[i].first is a Subsplit, with its bitset.first being the bitset representation of its clade Y
    if ( my_children[i].first.getYBitset().getFirstSetBit() == fsb || my_children[i].first.getZBitset().getFirstSetBit() == fsb )
    {
      if (u < my_children[i].second)
      {
        index = i;
        break;
      }
      u -= my_children[i].second;
    }
  }

  return my_children[index].first;
}

/* Draws a subsplit for subsplit S's clade Z*/
Subsplit SBNParameters::drawSubsplitForZ( const Subsplit &s ) const
{
  // Find all potential children of s
  std::vector<std::pair<Subsplit,double> > my_children = subsplit_cpds.at(s);

  // Find a distinguishing feature of clade Y in subsplit s
  // Since Y and Z are disjoint, we can use the first set bit in Z
  size_t fsb = s.getZBitset().getFirstSetBit();

  double u = GLOBAL_RNG->uniform01();
  size_t index;

  // This is like drawing a root split, but we must ensure we only draw from subsplits of Z
  for (size_t i=0; i<my_children.size(); ++i)
  {
    // This is a subsplit of Y if one of its splits has the same first set bit as Y
    // my_children[i].first is a Subsplit, with its bitset.first being the bitset representation of its clade Y
    if ( my_children[i].first.getYBitset().getFirstSetBit() == fsb || my_children[i].first.getZBitset().getFirstSetBit() == fsb )
    {
      if (u < my_children[i].second)
      {
        index = i;
        break;
      }
      u -= my_children[i].second;
    }
  }

  return my_children[index].first;
}

/* Checks the validity of an SBNParameters object.
   This requires checking that, for each candidate subsplit,
     1) The CPDs sum to 1 for both splitting Y and Z
     2) All candidate subsplits are compatible with their parents
     3) All subsplits are pairs of disjoint splits
   There are equivalent conditions on the root split
*/
bool SBNParameters::isValid(void) const
{
  if ( !isValidRootDistribution() )
  {
    return false;
  }

  std::pair<Subsplit,std::vector<std::pair<Subsplit,double> > > this_cpd;
  // Loop over parent subsplits
  BOOST_FOREACH(this_cpd, subsplit_cpds) {
    if ( !(isValidCPD(this_cpd.second, this_cpd.first)) )
    {
      return false;
    }
  }

  return true;
}

// Counts all subsplits in an unrooted tree (handles all virtual rooting)
void SBNParameters::countAllSubsplits(Tree& tree, std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts, std::map<Subsplit,double>& root_split_counts, std::map<Subsplit,double>& q, bool doSA)
{
  // Prep for tip to root pass
  std::string order = "postorder";
  tree.orderNodesForTraversal(order);
  const std::vector<TopologyNode*> &postorder_nodes = tree.getNodes();

  double one_over_n_branches = 1.0 / (2.0 * tree.getNumberOfNodes() - 3.0); // 1 over the number of branches in an unrooted tree

  // For storing subsplits
  std::vector<Subsplit> per_node_subsplit = std::vector<Subsplit>(tree.getNumberOfNodes(),Subsplit());
  // For storing sum of rooting probabilities from tip through this branch.
  std::vector<double> ttr = std::vector<double>(tree.getNumberOfNodes(),0.0);

  // Tip to root pass, here we do two things
  // 1) Get all nodes' subsplits (we will need these repeatedly)
  //      We do not make a subsplit for the root as the root is a trifurcation in an unrooted tree
  // 2) accumulate cumulative rooting probabilities and add to root split counters
  //      On each edge (the edge subtending the node we're visiting) we need the probability of rooting to the split this edge induces
  //      We do not loop over the root for now because it is a special case: we already account for rooting on all edges, but ttr[root] is ill-defined
  for (std::vector<TopologyNode*>::const_iterator it = postorder_nodes.begin(); it != (postorder_nodes.end()-1); ++it)
  {
    size_t index = (*it)->getIndex();

    if ( (*it)->isTip() )
    {
      // 1)
      per_node_subsplit[index] = (*it)->getSubsplit(taxa);

      // 2)
      Subsplit root_on_edge = per_node_subsplit[index].rootSplitFromClade();
      double this_root_q = doSA ? one_over_n_branches : q[root_on_edge];
      ttr[index] = this_root_q;
      incrementRootSplitCount(root_split_counts,root_on_edge,this_root_q);
    }
    else
    {
      std::vector<int> children = (*it)->getChildrenIndices();

      // 1)
      RbBitSet clade_1 = per_node_subsplit[children[0]].asCladeBitset();
      RbBitSet clade_2 = per_node_subsplit[children[1]].asCladeBitset();
      per_node_subsplit[index] = Subsplit(clade_1,clade_2);

      // 2)
      Subsplit root_on_edge = per_node_subsplit[index].rootSplitFromClade();
      double this_root_q = doSA ? one_over_n_branches : q[root_on_edge];
      ttr[index] = this_root_q + ttr[children[0]] + ttr[children[1]];
      incrementRootSplitCount(root_split_counts,root_on_edge,this_root_q);
    }
  }

  // Root to tip pass (this is where the fun starts)
  order = "preorder";
  tree.orderNodesForTraversal(order);
  const std::vector<TopologyNode*> &preorder_nodes = tree.getNodes();

  // Loop over edges of tree (exploit equivalency between an edge and the node that edge subtends)
  // The root has no edge so there is nothing to do for the root, so we skip it
  for (std::vector<TopologyNode*>::const_iterator it = preorder_nodes.begin()+1; it != preorder_nodes.end(); ++it)
  {
    // std::cout << ">>>working on a root/internal/tip node " << ((*it)->isRoot()) << "/" << ((*it)->isInternal()) << "/" << ((*it)->isTip()) << std::endl;
    // std::cout << ">The node's subsplit is " << per_node_subsplit[(*it)->getIndex()] << std::endl;
    // std::cout << ">The node's parent's subsplit is " << per_node_subsplit[(*it)->getParent().getIndex()] << std::endl;

    size_t index = (*it)->getIndex();

    // Edges descending from root need to be handled differently
    if ( (*it)->getParent().isRoot() )
    {
      // Get subsplits for other two descendants of root
      std::vector<int> root_children_indices = tree.getRoot().getChildrenIndices();

      std::vector<int> sibling_indices;

      for (size_t i=0; i<3; ++i)
      {
        if (index != root_children_indices[i])
        {
          sibling_indices.push_back(root_children_indices[i]);
        }
      }

      // Get all cases for virtual rooting of this edge (including current rooting)
      std::vector<std::pair<Subsplit,Subsplit> > cases = per_node_subsplit[index].doVirtualRootingRootParent(per_node_subsplit[sibling_indices[0]],per_node_subsplit[sibling_indices[1]],per_node_subsplit[index]);

      // Subsplit root_on_edge = per_node_subsplit[index].rootSplitFromClade();

      // Case 1
      double weight = ttr[root_children_indices[0]];
      incrementParentChildCount(parent_child_counts,cases[0],weight);
// std::cout << "did case 1" << std::endl;

      // Case 2
      weight = ttr[root_children_indices[0]];
      incrementParentChildCount(parent_child_counts,cases[1],weight);
// std::cout << "did case 2" << std::endl;

      // Case 3
      // weight = doSA ? one_over_n_branches : q[root_on_edge];
      weight = doSA ? one_over_n_branches : q[cases[5].first];
      incrementParentChildCount(parent_child_counts,cases[2],weight);
// std::cout << "did case 3" << std::endl;

      if ( !((*it)->isTip()) )
      {
        std::vector<int> children_indices = (*it)->getChildrenIndices();
        bool child_0_is_y = per_node_subsplit[index].isChildOfY(per_node_subsplit[children_indices[0]]);

        // Case 4
        weight = ttr[children_indices[child_0_is_y ? 0 : 1]];
        incrementParentChildCount(parent_child_counts,cases[3],weight);
// std::cout << "did case 4" << std::endl;

        // Case 5
        weight = ttr[children_indices[child_0_is_y ? 1 : 0]];
        incrementParentChildCount(parent_child_counts,cases[4],weight);
// std::cout << "did case 5" << std::endl;
      }

      // Case 6
      // weight = doSA ? one_over_n_branches : q[root_on_edge];
      weight = doSA ? one_over_n_branches : q[cases[5].first];
      incrementParentChildCount(parent_child_counts,cases[5],weight);
// std::cout << "did case 6" << std::endl;
    }
    else
    {
      // Define parent-child pair for current rooting (parent first, child second)
      std::pair<Subsplit,Subsplit> this_parent_child;
      this_parent_child.first = per_node_subsplit[(*it)->getParent().getIndex()];
      this_parent_child.second = per_node_subsplit[index];

      // Get all cases for virtual rooting of this edge (including current rooting)
      std::vector<std::pair<Subsplit,Subsplit> > cases = per_node_subsplit[index].doVirtualRootingNonRootParent(this_parent_child.first,this_parent_child.second);

      // Subsplit root_on_edge = per_node_subsplit[index].rootSplitFromClade();

      // Case 1
      double weight = 1 - ttr[(*it)->getParent().getIndex()];
      incrementParentChildCount(parent_child_counts,cases[0],weight);
// std::cout << "did case 1" << std::endl;

      // Case 2
      std::vector<int> my_parents_children = (*it)->getParent().getChildrenIndices();
      size_t sibling = 0;
      if (index == my_parents_children[0])
      {
        sibling = 1;
      }

      weight = ttr[sibling];
      incrementParentChildCount(parent_child_counts,cases[1],weight);
// std::cout << "did case 2" << std::endl;

      // Case 3
      // weight = doSA ? one_over_n_branches : q[root_on_edge];
      weight = doSA ? one_over_n_branches : q[cases[5].first];
      incrementParentChildCount(parent_child_counts,cases[2],weight);
// std::cout << "did case 3" << std::endl;

      if ( !((*it)->isTip()) )
      {
        std::vector<int> children_indices = (*it)->getChildrenIndices();
        bool child_0_is_y = per_node_subsplit[index].isChildOfY(per_node_subsplit[children_indices[0]]);

        // Case 4
        weight = ttr[children_indices[child_0_is_y ? 0 : 1]];
        incrementParentChildCount(parent_child_counts,cases[3],weight);
// std::cout << "did case 4" << std::endl;

        // Case 5
        weight = ttr[children_indices[child_0_is_y ? 1 : 0]];
        incrementParentChildCount(parent_child_counts,cases[4],weight);
// std::cout << "did case 5" << std::endl;
      }

      // Case 6
      // weight = doSA ? one_over_n_branches : q[root_on_edge];
      weight = doSA ? one_over_n_branches : q[cases[5].first];
      incrementParentChildCount(parent_child_counts,cases[5],weight);
// std::cout << "did case 6" << std::endl;

    }

  }

}

void SBNParameters::fitBranchLengthDistributions(std::vector<Tree> &trees )
{
  std::map<RbBitSet,std::vector<double> > branch_length_observations;

  // Loop over all trees
  for (size_t i=0; i<trees.size(); ++i)
  {
    double tree_length = 0.0;
    // Get branch lengths
    //unrooted trees have basal polytomies, so without fear we can take a node's clade, turn it into the split the edge below it represents, and add to that edge's observations
    const std::vector<TopologyNode*> tree_nodes = trees[i].getNodes();
    for (size_t i=0; i<tree_nodes.size(); ++i)
    {
      if ( (!tree_nodes[i]->isRoot()) )
      {
        Subsplit this_subsplit = tree_nodes[i]->getSubsplit(taxa);
        RbBitSet this_split = this_subsplit.asSplitBitset();

        (branch_length_observations[this_split]).push_back(tree_nodes[i]->getBranchLength());
        tree_length += tree_nodes[i]->getBranchLength();
      }
    }
    // Tree length
    RbBitSet whole_tree = RbBitSet(trees[i].getNumberOfTips(),true);
    (branch_length_observations[whole_tree]).push_back(tree_length);
  }

  if ( branch_length_approximation_method == "lognormalML" )
  {
    // Turn branch length observations into lognormal distributions
    std::pair<RbBitSet,std::vector<double> > clade_edge_observations;
    BOOST_FOREACH(clade_edge_observations, branch_length_observations) {
      std::vector<double> these_params = std::vector<double>(3,0.0);

      // std::cout << "Learning branch distribution for clade " << clade_edge_observations.first << ", observations are:" << std::endl;
      if (clade_edge_observations.second.size() > 2)
      {
        // TODO: if we are going to keep using lognormal via MLE, there are more efficient ways to get the logmean and logsd
        // Get mean/sd of log of observations
        double log_mean;
        for (size_t i=0; i<clade_edge_observations.second.size(); ++i)
        {
          // Branch lengths x such that c++ returns log(x) = -inf are possible, replace with smallest representable number instead
          double log_x = log(clade_edge_observations.second[i]);
          log_mean += log_x == RbConstants::Double::neginf ? RbConstants::Double::min : log_x;
        }
        log_mean /= clade_edge_observations.second.size();

        double log_sd;
        for (size_t i=0; i<clade_edge_observations.second.size(); ++i)
        {
          double log_x = log(clade_edge_observations.second[i]);
          log_sd += log_x == RbConstants::Double::neginf ? pow(RbConstants::Double::min - log_mean,2.0) : pow(log_x - log_mean,2.0);
        }
        log_sd /= clade_edge_observations.second.size();
        log_sd = sqrt(log_sd);

        // Approximate edge-length distribution using lognormal, use MLE parameters
        these_params[0] = log_mean;
        these_params[1] = log_sd;

        edge_length_distribution_parameters[clade_edge_observations.first] = these_params;

      }
      else
      {
        // Basically no information on edge length distribution
        // Approximate edge-length distribution using a lognormal resembling an exponential(10)
        std::vector<double> these_params = std::vector<double>(2,0.0);
        these_params[0] = -2.8;
        these_params[1] = 1.0;

        edge_length_distribution_parameters[clade_edge_observations.first] = these_params;
      }
    }
  }
  else if ( branch_length_approximation_method == "lognormalMOM" )
  {
    // Turn branch length observations into some distributions
    std::pair<RbBitSet,std::vector<double> > clade_edge_observations;

    BOOST_FOREACH(clade_edge_observations, branch_length_observations)
    {
      std::vector<double> these_params = std::vector<double>(2,0.0);
      if (clade_edge_observations.second.size() > 2)
      {
        std::vector<double> moments = get_moments(clade_edge_observations.second);
        these_params[0] = log(moments[0]/sqrt(1 + moments[1]/(pow(moments[0],2.0))));
        these_params[1] = sqrt(log(1 + moments[1]/(pow(moments[0],2.0))));
      }
      else
      {
        // Basically no information on edge length distribution
        // Approximate edge-length distribution using a lognormal resembling an exponential(10)
        these_params[0] = -2.8;
        these_params[1] = 1.0;
      }
      edge_length_distribution_parameters[clade_edge_observations.first] = these_params;
    }
  }
  else if ( branch_length_approximation_method == "generalizedGamma" )
  {
    // Turn branch length observations into some distributions
    std::pair<RbBitSet,std::vector<double> > clade_edge_observations;

    BOOST_FOREACH(clade_edge_observations, branch_length_observations)
    {
      std::vector<double> these_params = std::vector<double>(3,0.0);
      if (clade_edge_observations.second.size() > 100)
      {
        std::vector<double> gg_params = std::vector<double>(3,-1.0);
        gg_params = fit_gengam_gss(clade_edge_observations.second);

        // double gg_aic = 6 - fn_fixed_c_gengamma_lnl(gg_params[1],clade_edge_observations.second);
        // double gamma_aic = 4 - fn_fixed_c_gengamma_lnl(1.0,clade_edge_observations.second);
        //
        // // AIC-based comparison
        // if (gg_aic < gamma_aic)
        // {
          these_params = gg_params;
        // }
        // else
        // {
        //   std::vector<double> gamma_par = std::vector<double>(2,-1.0);
        //   gamma_par = fit_gamma_MOM(clade_edge_observations.second);
        //   these_params[0] = gamma_par[1];
        //   these_params[1] = 1.0;
        //   these_params[2] = gamma_par[0];
        // }

      }
      else if (clade_edge_observations.second.size() > 2)
      {
        std::vector<double> gamma_par = std::vector<double>(2,-1.0);
        gamma_par = fit_gamma_MOM(clade_edge_observations.second);
        these_params[0] = gamma_par[1];
        these_params[1] = 1.0;
        these_params[2] = gamma_par[0];
      }
      else
      {
        // Basically no information on edge length distribution
        // Approximate edge-length distribution using an exponential(10)
        these_params[0] = 10.0;
        these_params[1] = 1.0;
        these_params[2] = 1.0;
      }
      // std::cout << "fit gengamma to observations, params are a = " << these_params[0] << "; c = " << these_params[1] << "; l = " << these_params[2] << "; split was " << clade_edge_observations.first << std::endl;
      edge_length_distribution_parameters[clade_edge_observations.first] = these_params;
    }
  }
  else if ( branch_length_approximation_method == "gammaMOM" || branch_length_approximation_method == "compound" )
  {
    // Turn branch length observations into some distributions
    std::pair<RbBitSet,std::vector<double> > clade_edge_observations;
    BOOST_FOREACH(clade_edge_observations, branch_length_observations)
    {
      std::vector<double> these_params = std::vector<double>(2,0.0);

      if (clade_edge_observations.second.size() > 2)
      {
        these_params = fit_gamma_MOM(clade_edge_observations.second);
      }
      else
      {
        // Basically no information on edge length distribution
        // Approximate edge-length distribution using a lognormal resembling an exponential(10)
        these_params[0] = 1.0;
        these_params[1] = 10.0;
      }
      edge_length_distribution_parameters[clade_edge_observations.first] = these_params;
    }
  }

}

void SBNParameters::fitNodeTimeDistributions(std::vector<Tree> &trees )
{
  // First we get the data we want from our samples
  std::map<RbBitSet,std::vector<double> > node_time_observations;

  // Loop over all trees
  for (size_t i=0; i<trees.size(); ++i)
  {
    const std::vector<TopologyNode*> tree_nodes = trees[i].getNodes();

    // We need a mapping from the ith taxon in our bitsets to the corresponding tip node age
    // Our bitsets may not match bitsets available in tree, so we make our own mapping
    std::vector<double> tip_ages;

    for (size_t j=0; j<taxa.size(); ++j)
    {
      size_t index = trees[i].getTipIndex(taxa[j].getName());
      tip_ages.push_back(trees[i].getTipNode(index).getAge());
    }
// std::cout << "got tip ages" << std::endl;

    // Root is different, handle it first
    RbBitSet root_clade = RbBitSet(taxa.size(),true);

    double max_leaf_age = 0.0;

    for (size_t j=0; j<taxa.size(); ++j)
    {
      if ( tip_ages[j] > max_leaf_age)
      {
        max_leaf_age = tip_ages[j];
      }
    }

    double root_age = trees[i].getRoot().getAge();
// std::cout << "got root age" << std::endl;

    (node_time_observations[root_clade]).push_back(root_age - max_leaf_age);

    for (size_t n=0; n<tree_nodes.size(); ++n)
    {
      // Root already handled and nothing to do for tips
      if ( !(tree_nodes[n]->isTip()) && !(tree_nodes[n]->isRoot()) )
      {
// std::cout << "on node " << n << std::endl;
        Subsplit this_subsplit = tree_nodes[n]->getSubsplit(taxa);
        RbBitSet this_clade = this_subsplit.asCladeBitset();

        double max_descendant_leaf_age = 0.0;

        for (size_t j=0; j<taxa.size(); ++j)
        {
          if ( this_clade[j] && tip_ages[j] > max_descendant_leaf_age)
          {
            max_descendant_leaf_age = tip_ages[j];
          }
        }
// std::cout << "max_descendant_leaf_age = " << max_descendant_leaf_age << std::endl;
        double my_age_above_tips = tree_nodes[n]->getAge() - max_descendant_leaf_age;
// std::cout << "my_age_above_tips = " << my_age_above_tips << std::endl;
        double my_parents_age_above_tips = tree_nodes[n]->getParent().getAge() - max_descendant_leaf_age;
// std::cout << "my_parents_age_above_tips = " << my_parents_age_above_tips << std::endl;
        (node_time_observations[this_clade]).push_back(my_age_above_tips/my_parents_age_above_tips);
      }
    }
  }
// std::cout << "got all node ages" << std::endl;
  // Then we fit distributions to our observations
  if ( branch_length_approximation_method == "rootGammaNodePropKumaraswamy" )
  {
    // Turn root age observations into gamma distributions via method of moments
    RbBitSet root_clade = RbBitSet(taxa.size(),true);
    std::vector<double> root_age_observations = node_time_observations[root_clade];

    // Get mean/var of observations
    double mean;
    for (size_t i=0; i<root_age_observations.size(); ++i)
    {
      mean += root_age_observations[i];
    }
    mean /= root_age_observations.size();

    double var;
    for (size_t i=0; i<root_age_observations.size(); ++i)
    {
      var += pow(root_age_observations[i] - mean,2.0);
    }
    var /= root_age_observations.size();

    // Approximate edge-length distribution using gamma
    std::vector<double> root_params = std::vector<double>(2,0.0);

    root_params[1] = mean/var;
    root_params[0] = mean * root_params[1];

    edge_length_distribution_parameters[root_clade] = root_params;

    // Turn node age proportion observations into Kumaraswamy distributions
    std::pair<RbBitSet,std::vector<double> > node_proportion_observation;
    BOOST_FOREACH(node_proportion_observation, node_time_observations) {
      if ( node_proportion_observation.first != root_clade ) // We already did the root
      {
        if (node_proportion_observation.second.size() > 2)
        {
          // Approximate node-proportion distribution using Kumaraswamy
          edge_length_distribution_parameters[node_proportion_observation.first] = fit_kumar_agd(node_proportion_observation.second);

        }
        else
        {
          // Basically no information on node proportion distribution
          // Approximate node proportion distribution using a uniform
          std::vector<double> these_params = std::vector<double>(2,0.0);
          these_params[0] = 1.0;
          these_params[1] = 1.0;

          edge_length_distribution_parameters[node_proportion_observation.first] = these_params;
        }
      }
    }
  }

}
void SBNParameters::makeCPDs(std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts)
{

  // Put parent-child splits in correct format and place
  std::pair<std::pair<Subsplit,Subsplit>,double> this_parent_child;

  BOOST_FOREACH(this_parent_child, parent_child_counts) {
    Subsplit this_parent = this_parent_child.first.first;
    Subsplit this_child = this_parent_child.first.second;
    double this_prob = this_parent_child.second;

    if ( !(this_parent.isCompatible(this_child)) )
    {
      std::cout << "Invalid s|t" << std::endl;
      std::cout << "  s = " << this_child << std::endl;
      std::cout << "  t = " << this_parent << std::endl;
      throw(RbException("Found impossible parent-child subsplit pair in makeCPDs."));
    }

    std::pair<Subsplit,double> this_cpd;
    this_cpd.first = this_child;
    this_cpd.second = this_prob;

    (subsplit_cpds[this_parent]).push_back(this_cpd);
  }

  // Normalize CPDs
  std::pair<Subsplit,std::vector<std::pair<Subsplit,double> > > parent_cpd_pair;

  // Loop over parent subsplits
  BOOST_FOREACH(parent_cpd_pair, subsplit_cpds) {
    Subsplit parent = parent_cpd_pair.first; // The parent subsplit
    std::vector<std::pair<Subsplit,double> > my_children = parent_cpd_pair.second; // The children of this parent

    for (size_t i=0; i<my_children.size(); ++i)
    {
      if ( !(parent.isCompatible(my_children[i].first)) )
      {
        std::cout << "Found incompatible parent-child subsplit pair:" << parent << "->" << my_children[i].first << std::endl;
        throw(RbException("Found incompatible subsplit in makeCPDs."));
      }
    }

    normalizeCPDForSubsplit(my_children, parent);

  }

}

/*
  This function first puts the root split counts into a root splits probability map,
    then it normalizes them into a probability distribution
*/
void SBNParameters::makeRootSplits(std::map<Subsplit,double>& root_split_counts)
{
  // Put root splits in correct format and place
  std::pair<Subsplit,double> this_root;
  BOOST_FOREACH(this_root, root_split_counts) {
    root_splits.push_back(this_root);
  }

  // Normalize root splits
  double sum_root = 0.0;
  for (size_t i=0; i<root_splits.size(); ++i)
  {
    sum_root += root_splits[i].second;
  }

  for (size_t i=0; i<root_splits.size(); ++i)
  {
    root_splits[i].second /= sum_root;
  }

}

void SBNParameters::normalizeCPDForSubsplit(std::vector<std::pair<Subsplit,double> >& cpd, Subsplit& parent)
{

  double sum_y = 0.0; // sum of counts for child of parent subsplit's clade Y
  double sum_z = 0.0; // sum of counts for child of parent subsplit's clade Z

  // Find a distinguishing feature of clade Y in subsplit s
  // Since Y and Z are disjoint, we can use the first set bits in Y and Z
  size_t fsb_y = parent.getYBitset().getFirstSetBit();
  size_t fsb_z = parent.getZBitset().getFirstSetBit();

  // In unrooted counting, we lose dummy subsplits, we add them in here
  bool y_is_tip = parent.getYBitset().getNumberSetBits() == 1 ? true : false;
  bool z_is_tip = parent.getZBitset().getNumberSetBits() == 1 ? true : false;

  size_t n_children_of_y = 0;
  size_t n_children_of_z = 0;

  for (size_t i=0; i<cpd.size(); ++i) // Loop over the children of this parent, get sum for normalizing
  {
    if ( !(parent.isCompatible(cpd[i].first)) )
    {
      std::cout << "Found incompatible parent-child subsplit pair." << std::endl;
      std::cout << "  parent: " << parent << std::endl;
      std::cout << "  child:  " << cpd[i].first << std::endl;
      throw(RbException("Found incompatible subsplit in normalizeCPDForSubsplit."));
    }

    // This is a subsplit of parent's clade Y if one of its splits has the same first set bit as Y
    // cpd[i].first is a Subsplit, with its bitset.first being the bitset representation of its clade Y
    if ( cpd[i].first.getYBitset().getFirstSetBit() == fsb_y || cpd[i].first.getZBitset().getFirstSetBit() == fsb_y )
    {
      sum_y +=  cpd[i].second;
      ++n_children_of_y;
    }
    else if ( cpd[i].first.getYBitset().getFirstSetBit() == fsb_z || cpd[i].first.getZBitset().getFirstSetBit() == fsb_z )
    {
      sum_z +=  cpd[i].second;
      ++n_children_of_z;
    }
    else {
      throw(RbException("Found incompatible subsplit when learning SBN."));
    }
  }

  // If y or z are tips, we will now insert probability-1 dummy splits for these tips
  if ( n_children_of_y == 0 && y_is_tip )
  {
    RbBitSet tip = parent.getYBitset();
    Subsplit dummy_subsplit = Subsplit(tip,tip);

    std::pair<Subsplit,double> dummy_cpd;
    dummy_cpd.first = dummy_subsplit;
    dummy_cpd.second = 1.0;

    (subsplit_cpds[parent]).push_back(dummy_cpd);
    sum_y = 1.0;
  }

  if ( n_children_of_z == 0 && z_is_tip )
  {
    RbBitSet tip = parent.getZBitset();
    Subsplit dummy_subsplit = Subsplit(tip,tip);

    std::pair<Subsplit,double> dummy_cpd;
    dummy_cpd.first = dummy_subsplit;
    dummy_cpd.second = 1.0;

    (subsplit_cpds[parent]).push_back(dummy_cpd);
    sum_z = 1.0;
  }

  for (size_t i=0; i<cpd.size(); ++i) // Loop over the children of this parent, normalize
  {
    // This is a subsplit of X's clade Y if one of its splits has the same first set bit as Y
    // cpd[i].first is a Subsplit, with its bitset.first being the bitset representation of its clade Y
    if ( cpd[i].first.getYBitset().getFirstSetBit() == fsb_y || cpd[i].first.getZBitset().getFirstSetBit() == fsb_y )
    {
      (subsplit_cpds[parent][i]).second /= sum_y;
    }
    else
    {
      (subsplit_cpds[parent][i]).second /= sum_z;
    }
  }

}

bool SBNParameters::isValidCPD(std::vector<std::pair<Subsplit,double> >& cpd, Subsplit& parent) const
{
  double sum_y = 0.0;
  double sum_z = 0.0;

  size_t fsb_y = parent.getYBitset().getFirstSetBit();
  size_t fsb_z = parent.getZBitset().getFirstSetBit();

  bool y_is_tip = parent.getYBitset().getNumberSetBits() == 1 ? true : false;
  bool z_is_tip = parent.getZBitset().getNumberSetBits() == 1 ? true : false;

  for (size_t i=0; i<cpd.size(); ++i)
  {
    Subsplit child = cpd[i].first;
    size_t child_fsb = child.asCladeBitset().getFirstSetBit();
    RbBitSet child_y = child.getYBitset();
    RbBitSet child_z = child.getZBitset();

    // Make sure child is compatible with parent
    if ( !(parent.isCompatible(child)) )
    {
      std::cout << "Found incompatible parent-child subsplit pair." << std::endl;
      std::cout << "  parent: " << parent << std::endl;
      std::cout << "  child:  " << child << std::endl;
      return false;
    }

    // Make sure children's splits are disjoint
    if ( !child.splitsAreDisjoint() )
    {
      std::cout << "Found impossible subsplit: " << child << std::endl;
      return false;
    }

    // Determine if the child is of parent's split Y or Z, add probability to appropriate sum
    if ( child_fsb == fsb_y )
    {
      sum_y += cpd[i].second;
    }
    else
    {
      sum_z += cpd[i].second;
    }
  }

  // Make sure that for all non-trivial subsplits, the sum of probabilities of children is 1 (trivial subsplits are tips and don't need CPDs)
  double tol = 0.0001;
  // if ( (fabs(sum_y - 1.0) > tol && !y_is_tip) || (fabs(sum_z - 1.0) > tol && !z_is_tip) )
  if ( (fabs(sum_y - 1.0) > tol) || (fabs(sum_z - 1.0) > tol) )
  {
    std::cout << "Unnormalized or improperly normalized CPD for parent subsplit " << parent << std::endl;
    std::cout << "Sum of CPDs for descendants of Y is " << std::fixed << std::setprecision(10) << sum_y << ". Sum of CPDs for descendants of Z is " << sum_z << std::endl;
    return false;
  }

  return true;
}

bool SBNParameters::isValidRootDistribution(void) const
{
  double sum_root = 0.0;
  for (size_t i=0; i<root_splits.size(); ++i)
  {
    sum_root += root_splits[i].second;
    RbBitSet root_y = root_splits[i].first.getYBitset();
    RbBitSet root_z = root_splits[i].first.getZBitset();

    // Check that in sum_root they have all tips
    if ( !(root_y.getNumberSetBits() + root_z.getNumberSetBits() == root_y.size()) )
    {
      return false;
    }
    // Check they're disjoint
    if ( !root_splits[i].first.splitsAreDisjoint() )
    {
      std::cout << "Found impossible root split: " << root_splits[i].first << std::endl;
        return false;
    }
    // Check they're not empty
    if ( root_y.getNumberSetBits() == 0 || root_z.getNumberSetBits() == 0 )
    {
      std::cout << "Found impossible root split: " << root_splits[i].first << ". Pr(root_split) = " << root_splits[i].second << std::endl;
        return false;
    }
  }

  double tol = 0.0001;
  if ( fabs(sum_root - 1.0) > tol )
  {
    std::cout << "Root splits are unnormalized or improperly normalized, sum of root probabilities is " << sum_root << std::endl;
    return false;
  }

  return true;
}

void SBNParameters::learnTimeCalibratedSBN( std::vector<Tree>& trees )
{
  time_calibrated = true;

  if ( !(branch_length_approximation_method == "rootGammaNodePropKumaraswamy") )
  {
    throw(RbException("Invalid branch length/node height approximation method when initializing SBN object."));
  }

  // For counting subsplits, we could use integers but unrooted trees get fractional counts, so we'll be consistent
  std::map<Subsplit,double> root_split_counts;
  std::map<std::pair<Subsplit,Subsplit>,double> parent_child_counts;

  // The weight to assign when counting subsplits, for rooted trees the weight is 1
  double weight = 1.0;

  // Loop over all trees
  // for each, get all root splits and subsplit parent-child relationships
  // then consolidate into our master list
  for (size_t i=0; i<trees.size(); ++i)
  {
    addTreeToAllRootSplitCounts(root_split_counts, trees[i], weight);
    addTreeToAllParentChildCounts(parent_child_counts, trees[i], weight);
// std::cout << "got splits from tree " << i << std::endl;
  }

  // Turn root split counts into a distribution on the root split
  makeRootSplits(root_split_counts);
// std::cout << "made root splits" << std::endl;

  // Turn parent-child subsplit counts into CPDs
  makeCPDs(parent_child_counts);
// std::cout << "made CPDs" << std::endl;

  if ( !isValid() )
  {
    throw(RbException("learnTimeCalibratedSBN produced an invalid SBNParameters object."));
  }

  fitNodeTimeDistributions(trees);
// std::cout << "fit branch lengths" << std::endl;

}

void SBNParameters::learnUnconstrainedSBNSA( std::vector<Tree> &trees )
{
  time_calibrated = false;

  if ( !(branch_length_approximation_method == "generalizedGamma" || branch_length_approximation_method == "compound" || branch_length_approximation_method == "gammaMOM" || branch_length_approximation_method == "lognormalML" || branch_length_approximation_method == "lognormalMOM") )
  {
    throw(RbException("Invalid branch length/node height approximation method when initializing SBN object."));
  }

  // std::cout << "hello from learnUnconstrainedSBNSA, there are this many trees" << trees.size() << std::endl;

  // To store counts
  std::map<Subsplit,double> root_split_counts;
  std::map<std::pair<Subsplit,Subsplit>,double> parent_child_counts;

  // This can stay empty, we don't need to specify q() if we override with doSA=TRUE
  std::map<Subsplit,double> q;

  // Run counting
  for (size_t i=0; i<trees.size(); ++i)
  {
// std::cout << trees[i] << std::endl;
    countAllSubsplits(trees[i], parent_child_counts, root_split_counts, q, true);
  }
// std::cout << "counted all subsplits in all trees" << std::endl;

  // Turn root split counts into a distribution on the root split
  makeRootSplits(root_split_counts);

  // Turn parent-child subsplit counts into CPDs
  makeCPDs(parent_child_counts);

  // Handle branch lengths
  fitBranchLengthDistributions(trees);

  if ( !isValid() )
  {
    throw(RbException("learnUnconstrainedSBNSA produced an invalid SBNParameters object."));
  }

}

void SBNParameters::learnUnconstrainedSBNEM( std::vector<Tree> &trees, double &alpha )
{

  // Initialize parameters to SA values
  // This will handle the branch lengths for us, we don't need to touch them again
  learnUnconstrainedSBNSA(trees);

  // // run EM, start on E-step because we have parameters from M-step
  // bool terminate = false;
  // while ( !terminate )
  // {
  //   // E-step, compute q and m
  //
  //   // M-step, compute p
  //
  // }
}

// Takes a tree in with weight (can be for multiple trees or for q(root)), keeps rooting intact, adds all parent-child subsplits found in this tree to master list of counts
void SBNParameters::addTreeToAllParentChildCounts(std::map<std::pair<Subsplit,Subsplit>,double> &parent_child_counts, Tree& tree, double &weight)
{
  std::vector<std::pair<Subsplit,Subsplit> > these_parent_child_subsplits = tree.getAllSubsplitParentChildPairs(taxa);
  for (size_t j=0; j<these_parent_child_subsplits.size(); ++j)
  {
    if ( parent_child_counts.count(these_parent_child_subsplits[j]) == 0 )
    {
      parent_child_counts[these_parent_child_subsplits[j]] = weight;
    }
    else
    {
      parent_child_counts[these_parent_child_subsplits[j]] += weight;
    }
  }
}

// Takes a tree in with weight (can be for multiple trees or for q(root)), keeps rooting intact, adds the root split to master list of counts
void SBNParameters::addTreeToAllRootSplitCounts(std::map<Subsplit,double>& root_split_counts, Tree& tree, double &weight)
{
  Subsplit this_root_split;
  this_root_split = tree.getRootSubsplit(taxa);
  if ( root_split_counts.count(this_root_split) == 0 )
  {
    root_split_counts[this_root_split] = weight;
  }
  else
  {
    root_split_counts[this_root_split] += weight;
  }
}


void SBNParameters::incrementParentChildCount(std::map<std::pair<Subsplit,Subsplit>,double> &parent_child_counts, std::pair<Subsplit,Subsplit> &this_parent_child, double &weight)
{
// std::cout << "incrementing ParentChildCount by " << weight << " for parent-child" << this_parent_child.first << " - " << this_parent_child.second << std::endl;
  if ( parent_child_counts.count(this_parent_child) == 0 )
  {
    parent_child_counts[this_parent_child] = weight;
  }
  else
  {
    parent_child_counts[this_parent_child] += weight;
  }
}

void SBNParameters::incrementRootSplitCount(std::map<Subsplit,double>& root_split_counts, Subsplit &this_root_split, double &weight)
{
// std::cout << "incrementing RootSplitCount by " << weight << " for root split" << this_root_split << std::endl;
  if ( root_split_counts.count(this_root_split) == 0 )
  {
    root_split_counts[this_root_split] = weight;
  }
  else
  {
    root_split_counts[this_root_split] += weight;
  }
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const SBNParameters& x) {

    // std::streamsize previousPrecision = o.precision();
    // std::ios_base::fmtflags previousFlags = o.flags();

    std::vector<Taxon> taxa = x.getTaxa();
    o << "SBN on taxon vector [ ";
    // o << std::fixed;
    // o << std::setprecision(4);

    o << taxa[0].getName();

    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=1; i < x.getNumTaxa(); i++)
    {
      o << ", ";
      o << taxa[i].getName();
    }

    o << " ]";
    // o.setf(previousFlags);
    // o.precision(previousPrecision);

    return o;
}
