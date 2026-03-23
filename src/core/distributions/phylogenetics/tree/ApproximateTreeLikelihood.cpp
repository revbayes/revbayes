#include "ApproximateTreeLikelihood.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>

#include "DistributionMultivariateNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
    topology_match_checked( false )
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
    // Check that the topologies are identical

    // make the time tree unrooted
    const Tree &time_tree_copy = time_tree->getValue();
    Tree* time_tree_unrooted = time_tree_copy.clone();
    time_tree_unrooted->unroot();
        
    std::string newick_time_tree = time_tree_unrooted->getPlainNewickRepresentation();
    std::string newick_branch_length_tree = value->getPlainNewickRepresentation();

    return newick_time_tree == newick_branch_length_tree;
}


double ApproximateTreeLikelihood::computeLnProbability( void )
{

    // only check once if the topologies are the same
    if ( topology_match_checked )
    {
        bool match = checkTopologyMatch();
        if ( match == false )
        {
            return RbConstants::Double::neginf;
        }
        topology_match_checked = true;
    }
    
    // get all our variables to compute the approximate likelihood
    std::vector<double> branch_lengths = computeBranchLengths();
    size_t num_branches = gradients->size();

    double ln_prob = 0.0;
    for (size_t i=0; i<num_branches; ++i)
    {
        // std::cerr << "The " << i << "-th branch length is " << mle_branch_lengths[i] << " and the " << i << "-th gradient is " << (*gradients)[i] << std::endl;
        ln_prob += mle_branch_lengths[i] * (*gradients)[i];
        ln_prob += mle_branch_lengths[i] * mle_branch_lengths[i] * (*hessian)[i][i] / 2.0;
        for (size_t j=0; i<j; ++j)
        {
            ln_prob += mle_branch_lengths[i] * mle_branch_lengths[j] * (*hessian)[i][j];
        }
    }
    
    return ln_prob;
}


std::vector<double> ApproximateTreeLikelihood::computeBranchLengths( void )
{
    
    // get time tree nodes
    const std::vector<TopologyNode*> &time_tree_nodes = time_tree->getValue().getNodes();
    size_t num_tips = time_tree->getValue().getNumberOfTips();
    
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
            
            // get the matching index to be sorted correctly
            RbBitSet split = RbBitSet(num_tips);
            the_node->getTaxa(split);
            
            size_t index = split_to_index[split];
            
            // we use the plus because the two root branches need to be added
            // otherwise the plus actually never does anything
            branch_lengths[index] += new_branch_length;
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
    
    size_t num_tips = value->getNumberOfTips();
    split_to_index.clear();
    
    // loop through tree nodes and draw new rates to get new branch lengths
    for (size_t i=0; i<bl_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = bl_tree_nodes[i];
        // std::cerr << "The branch length of the " << i << "-th node is " << the_node->getBranchLength() << std::endl;
                
        RbBitSet split = RbBitSet(num_tips);
        the_node->getTaxa(split);
        
        RbBitSet split_rev = split;
        split_rev.flip();
        
        split_to_index[split] = i;
        split_to_index[split_rev] = i;

        // get the branch time
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
    
    // simulate displacement around the MLE
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
    split_to_index.clear();
    
    // loop through tree nodes and draw new rates to get new branch lengths
    for (size_t i=0; i<bl_tree_nodes.size(); ++i)
    {
        TopologyNode* the_node = bl_tree_nodes[i];
        
        RbBitSet split = RbBitSet(num_tips);
        the_node->getTaxa(split);
        
        RbBitSet split_rev = split;
        split_rev.flip();
        
        split_to_index[split] = i;
        split_to_index[split_rev] = i;

        // get the branch time
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
