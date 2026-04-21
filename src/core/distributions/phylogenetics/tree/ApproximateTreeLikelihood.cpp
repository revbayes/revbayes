#include "ApproximateTreeLikelihood.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>

#include "DistributionMultivariateNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbBitSet.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TreeChangeEventMessage.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

ApproximateTreeLikelihood::ApproximateTreeLikelihood(const TypedDagNode<Tree>* tt, TypedDagNode< RbVector<double> > *br, const RbVector<double> *gr, const MatrixReal *h, TRANSFORMATION tr) : TypedDistribution<Tree>( new Tree() ),
    time_tree( tt ),
    branch_rates( br ),
    gradients( gr ),
    hessian( h ),
    transform( tr ),
    topology_match_checked( false ),
    mle_preorder( false )
{

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( branch_rates );
    this->addParameter( time_tree );

    simulateTree();

}


ApproximateTreeLikelihood::~ApproximateTreeLikelihood()
{
    
    // the tree will be deleted automatically by the base class

}


ApproximateTreeLikelihood* ApproximateTreeLikelihood::clone( void ) const
{

    return new ApproximateTreeLikelihood( *this );
}


bool ApproximateTreeLikelihood::checkTopologyMatch( void ) const
{

    return time_tree->getValue().hasSameTopology( *value );
}


/**
 * Fill internal_to_mle_index so thatcomputeBranchLengths() can add each proposed branch length into branch_lengths[j]
 * with j = internal_to_mle_index[i], matching the indexing of mle_branch_lengths[j], gradients[j], and Hessian entries
 * that use j. This mapping step is necessary because the order of time_tree.getNodes() is not the same as the order
 * used to build mle_branch_lengths. Matching is done by tip bipartition (split): each non-root node defines the split
 * across the branch above it.
 *
 * Steps:
 *   (1) split_to_mle: bipartition -> index j in mle_branch_lengths.
 *   (2) For each non-root time-tree node i in getNodes() order, internal_to_mle_index[i] = split_to_mle[split].
 *
 * Both orientations of each split are inserted because getTaxa() may return either side of the bipartition.
 */
void ApproximateTreeLikelihood::initializeMapping( void )
{
    const Tree &tt = time_tree->getValue();
    const std::vector<TopologyNode*> &time_tree_nodes = tt.getNodes();
    const size_t num_tips = tt.getNumberOfTips();

    // We clear the mapping vector and re-initialize it with one slot per time-tree node index (same order as getNodes()).
    // We fill it with max() values to make a clear distinction between unset (unused) and set entries: this value cannot
    // be a valid branch index. The root gets to keep this value, since there is no branch above the root to map.
    internal_to_mle_index.clear();
    internal_to_mle_index.resize(time_tree_nodes.size(), std::numeric_limits<size_t>::max());

    std::map<RbBitSet, size_t> split_to_mle;

    // Step 1: find j such that j is the entry in `mle_branch_lengths` for the branch that carries the split defined by one
    // branch of `value`. Build `split_to_mle` to store this mapping.
    if ( mle_preorder == true )
    {
        // Matches setValue(): same preorder layout after reorderPreorder. Node index k runs 1, ..., N-1 (N = total nodes);
        // k-1 indexes mle_branch_lengths. k == 0 is the root and is omitted because only branches with a parent edge
        // appear in mle_branch_lengths. Therefore, the vector has length N-1 and indices 0, ..., N-2.
        Tree *v_copy = value->clone();
        TreeUtilities::reorderPreorder( *v_copy );
        for (size_t k = 1; k < v_copy->getNumberOfNodes(); ++k)
        {
            RbBitSet split = RbBitSet(num_tips);
            v_copy->getNode(k).getTaxa(split);
            RbBitSet split_rev = split;
            split_rev.flip();
            split_to_mle[split] = k - 1;
            split_to_mle[split_rev] = k - 1;
        }
        delete v_copy;
    }
    else
    {
        // Matches simulateTree(): MLE slots follow value->getNodes() order after unroot().
        const std::vector<TopologyNode*> &bl_nodes = value->getNodes();
        for (size_t i = 0; i < bl_nodes.size(); ++i)
        {
            RbBitSet split = RbBitSet(num_tips);
            bl_nodes[i]->getTaxa(split);
            RbBitSet split_rev = split;
            split_rev.flip();
            split_to_mle[split] = i;
            split_to_mle[split_rev] = i;
        }
    }

    // Step (2): for each non-root time-tree node i, record which mle_branch_lengths index j matches its split.
    for (size_t i = 0; i < time_tree_nodes.size(); ++i)
    {
        if ( time_tree_nodes[i]->isRoot() == true )
        {
            continue;
        }

        RbBitSet split = RbBitSet(num_tips);
        time_tree_nodes[i]->getTaxa(split);
        const auto it = split_to_mle.find(split);
        if ( it == split_to_mle.end() )
        {
            throw RbException("Could not match an internal branch to a corresponding branch in the MLE tree.");
        }
        internal_to_mle_index[i] = it->second;
    }
}


double ApproximateTreeLikelihood::computeLnProbability( void )
{

    // only check once if the topologies are the same
    if ( topology_match_checked == false )
    {
        bool match = checkTopologyMatch();
        if ( match == false )
        {
            return RbConstants::Double::neginf;
        }
        topology_match_checked = true;
        initializeMapping();
    }
    
    // get the current branch lengths from the MCMC state
    std::vector<double> branch_lengths = computeBranchLengths();
    size_t num_branches = gradients->size();

    // second-order Taylor expansion of ln L around the MLE branch lengths:
    // ln L(b) ≈ const + g'(b - b_mle) + 1/2 (b - b_mle)' H (b - b_mle)
    std::vector<double> delta(num_branches);
    double ln_prob = 0.0;

    for (size_t i = 0; i < num_branches; ++i)
    {
        const RbVector<double> &hessian_i = (*hessian)[i];
        delta[i] = branch_lengths[i] - mle_branch_lengths[i];

        ln_prob += delta[i] * (*gradients)[i];
        ln_prob += delta[i] * delta[i] * hessian_i[i] / 2.0;
        
        for (size_t j = 0; j < i; ++j)
        {
            ln_prob += delta[i] * delta[j] * hessian_i[j];
        }
    }
    
    return ln_prob;
}


std::vector<double> ApproximateTreeLikelihood::computeBranchLengths( void )
{
    
    // get time tree nodes
    const std::vector<TopologyNode*> &time_tree_nodes = time_tree->getValue().getNodes();
    
    // initialize the branch length vector: need one fewer nodes than there are in the (rooted) time tree
    std::vector<double> branch_lengths(time_tree_nodes.size()-1, 0.0);

    // loop through time tree nodes and draw new rates to get new branch lengths
    for (size_t i=0; i<time_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = time_tree_nodes[i];
        
        if ( the_node->isRoot() == false )
        {
            // get the branch time
            double branch_time = the_node->getBranchLength();

            // get branch rate
            double new_branch_rate = branch_rates->getValue()[i];
            
            // get new branch length
            double new_branch_length = new_branch_rate * branch_time;

            // find the index of the corresponding branch in the MLE tree
            const size_t mle_index = internal_to_mle_index[i];

            // we use the plus because the two root branches of the MLE tree are merged into one in the time tree
            // otherwise the plus actually never does anything
            branch_lengths[mle_index] += new_branch_length;
        }
    }

    return branch_lengths;
}


void ApproximateTreeLikelihood::transformBranchLengths( std::vector<double>& branch_lengths )
{
    size_t nb = branch_lengths.size();
    
    for (size_t i=0; i<nb; ++i)
    {
        if ( transform == TRANSFORMATION::LOG )
        {
            branch_lengths[i] = log(branch_lengths[i]);
        }
        else if ( transform == TRANSFORMATION::SQRT )
        {
            branch_lengths[i] = sqrt(branch_lengths[i]);
        }
        else if ( transform == TRANSFORMATION::ARCSIN )
        {
            branch_lengths[i] = asin(branch_lengths[i]);
        }
    }
}


void ApproximateTreeLikelihood::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{

    if ( m == TreeChangeEventMessage::TOPOLOGY )
    {
        throw RbException("Don't do that!!");
    }

}


void ApproximateTreeLikelihood::keepSpecialization(const DagNode* affecter)
{
    
}


void ApproximateTreeLikelihood::redrawValue( void )
{
    simulateTree();
}


void ApproximateTreeLikelihood::setValue(RevBayesCore::Tree *v, bool force)
{
    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );

    // Sebastian: check if anything special needs to be done
    
    // copy the branch-length tree, reorder the copy by preorder traversal, and get its nodes
    Tree* v_copy = v->clone();
    TreeUtilities::reorderPreorder( *v_copy );
    std::vector<TopologyNode*> bl_tree_nodes;
    
    for (size_t i = 1; i < v_copy->getNumberOfNodes(); ++i)  // skip root (i=0)
    {
        bl_tree_nodes.push_back( &v_copy->getNode(i) );
    }
    
    mle_preorder = true;
    mle_branch_lengths.resize( bl_tree_nodes.size() );

    // loop through the branch-length tree nodes and store the new branch lengths
    for (size_t i=0; i<bl_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = bl_tree_nodes[i];
        // std::cerr << "The branch length of the " << i << "-th node is " << the_node->getBranchLength() << std::endl;
        mle_branch_lengths[i] = the_node->getBranchLength();
    }
    
    transformBranchLengths( mle_branch_lengths );
    topology_match_checked = false;
    delete v_copy;
}


void ApproximateTreeLikelihood::simulateTree( void )
{
    delete value;

    // make our own copy of this time tree
    // our strategy is to loop through all nodes of the time tree, update
    // the branch times with new branch lengths, and then finally unroot
    // the tree to create our new branch length tree
    const Tree &time_tree_copy = time_tree->getValue();
    value = time_tree_copy.clone();
    
    // we need to change the tree setting so that the branch lengths are used and not the ages.
    // internally, our tree nodes store both an age and branch length variable, but only one should be used.
    value->getRoot().setUseAges( false, true );

    // get time tree nodes
    const std::vector<TopologyNode*> &time_tree_nodes = time_tree->getValue().getNodes();
    size_t num_tips = time_tree->getValue().getNumberOfTips();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;

    // With approximate likelihood, we assume that the branch-length trees we are in effect sampling (by sampling
    // time trees and branch rates) are displaced around the MLE tree by a multivariate normal distribution. Here,
    // we assume by symmetry that we can simulate "MLE" trees using multivariate normal displacements around the
    // current proposed tree.
    std::vector<double> displacements = RbStatistics::MultivariateNormal::rvCovariance(static_cast<const std::vector<double>&>(*gradients), *hessian, *rng);

    // loop through time tree nodes and draw new rates to get new branch lengths
    for (size_t i=0; i<time_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = time_tree_nodes[i];

        if ( the_node->isRoot() == false )
        {
            // get the branch time
            double branch_time = the_node->getBranchLength();

            // get branch rate
            double new_branch_rate = branch_rates->getValue()[i];

            // get new branch length
            // SH: maybe add code here to transform branches (log, arcsin, sqrt)
            // David: Maybe make nicer
            double d = 0.0;
            if ( i < displacements.size() )
            {
                d = displacements[i];
            }
            double new_branch_length = new_branch_rate * branch_time + d;

            new_branch_length = fmax(0,new_branch_length);

            // set new branch length
            the_node->setBranchLength( new_branch_length );
        }
    }

    // now unroot to get the final branch length tree
    value->unroot();
    
    const std::vector<TopologyNode*> &bl_tree_nodes = value->getNodes();
    mle_branch_lengths = std::vector<double>(bl_tree_nodes.size(), 0.0);
    mle_preorder = false;

    // loop through the branch-length tree nodes and store the new branch lengths
    for (size_t i=0; i<bl_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = bl_tree_nodes[i];
        mle_branch_lengths[i] = the_node->getBranchLength();
    }
    
    transformBranchLengths( mle_branch_lengths );
    
    topology_match_checked = false;

}


void ApproximateTreeLikelihood::restoreSpecialization(const DagNode *restorer)
{
    
}


/** Swap a parameter of the distribution */
void ApproximateTreeLikelihood::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    if (oldP == branch_rates )
    {
        branch_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == time_tree )
    {
        time_tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }

}


void ApproximateTreeLikelihood::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    
}
