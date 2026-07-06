#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "DistributionMultivariateNormal.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "MatrixReal.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloMultivariateBrownianProcessREML::PhyloMultivariateBrownianProcessREML(const TypedDagNode<Tree> *t, const TypedDagNode<MatrixReal> *c, size_t ns) :
    AbstractPhyloBrownianProcess( t, ns ),
    node_likelihoods(this->num_nodes),
    independent_contrasts( std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0.0) ) ),
    independent_contrasts_sds( std::vector<double>(this->num_nodes, 0.0) ),
    rate_matrix( c ),
    active_matrix(0),
    precision_matrices( std::vector<MatrixReal>( 2, MatrixReal(num_sites) ) )
{
    
    // add the parameters to our set
    this->addParameter( rate_matrix );
    
    // make sure the rate matrix is inverted using Cholesky decomposition
    rate_matrix->getValue().setCholesky(true);
    
    // compute the inverse variance-covariance matrix (the precision matrix)
    precision_matrices[0] = rate_matrix->getValue().computeInverse();
    precision_matrices[0].setCholesky(true);
    precision_matrices[1] = rate_matrix->getValue().computeInverse();
    precision_matrices[1].setCholesky(true);
    
    // We don't want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );

    // now we need to reset the values
    this->redrawValue();
    
    // we need to reset the contrasts
    resetValue();

}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloMultivariateBrownianProcessREML::~PhyloMultivariateBrownianProcessREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloMultivariateBrownianProcessREML* PhyloMultivariateBrownianProcessREML::clone( void ) const
{
    
    return new PhyloMultivariateBrownianProcessREML( *this );
}


double PhyloMultivariateBrownianProcessREML::computeLnProbability( void )
{
    
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        resetValue();
    }
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( not this->node_likelihoods.is_valid(rootIndex) )
    {
        recursiveComputeLnProbability( root, rootIndex );
    }

    // return the likelihood at the root
    this->ln_prob = this->node_likelihoods[rootIndex].partial_likelihood;
        
    
    return this->ln_prob;
    
}



void PhyloMultivariateBrownianProcessREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}



std::vector<std::vector<double> > PhyloMultivariateBrownianProcessREML::getContrasts()
{

    // make sure the necessary quantities are clean
    computeLnProbability();
    
    // clean out the contrasts vector
    independent_contrasts.clear();
    independent_contrasts = std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0.0) );

    independent_contrasts_sds.clear();
    independent_contrasts_sds = std::vector<double>(this->num_nodes, 0.0);

    // compute the constrasts by recursively calling the contrast calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();

    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();

    // compute the contrast at the root, and recursively compute the remaining contrasts up the tree
    recursiveComputeContrasts( root, rootIndex );
    
    // remove the first n standardized contrasts, where n is the number of tips
    size_t num_tips = this->tau->getValue().getNumberOfTips();
    
    independent_contrasts.erase(independent_contrasts.begin(), independent_contrasts.begin() + num_tips);
    independent_contrasts_sds.erase(independent_contrasts_sds.begin(), independent_contrasts_sds.begin() + num_tips);

    return independent_contrasts;
    
}


void PhyloMultivariateBrownianProcessREML::keepSpecialization( const DagNode* affecter )
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.keep();
    }
}


void PhyloMultivariateBrownianProcessREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{

    // check for recomputation
    if ( node.isTip() == false && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihood = 0.0;
        node_cache.contrasts.assign(this->num_sites, 0.0);
        node_cache.contrast_uncertainty = 0.0;

        double              &p_node  = node_cache.partial_likelihood;
        std::vector<double> &mu_node = node_cache.contrasts;

        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        
        for (size_t j = 1; j < num_children; ++j)
        {
        
            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left = &node.getChild(0);
                left_index = left->getIndex();
                recursiveComputeLnProbability( *left, left_index );
            }

            const NodeCache *left_cache = (j == 1 ? &this->node_likelihoods[left_index] : &node_cache);
            
            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeLnProbability( right, right_index );

            const NodeCache &right_cache = this->node_likelihoods[right_index];
            const double &p_left  = left_cache->partial_likelihood;
            const double &p_right = right_cache.partial_likelihood;

            // get the per node and site contrasts
            const std::vector<double> &mu_left  = left_cache->contrasts;
            const std::vector<double> &mu_right = right_cache.contrasts;

            // get the propagated uncertainties
            double delta_left  = left_cache->contrast_uncertainty;
            double delta_right = right_cache.contrast_uncertainty;

            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());

            // add the propagated uncertainty to the branch lengths
            double t_left  = v_left  + delta_left;
            double t_right = v_right + delta_right;

            // set delta_node = (t_l*t_r)/(t_l+t_r);
            node_cache.contrast_uncertainty = (t_left * t_right) / (t_left + t_right);

            double branch_length = t_left + t_right;
            
            std::vector<double> these_contrasts(num_sites);
            std::vector<double> means(num_sites, 0.0);
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                // compute the contrasts for this site and node
                these_contrasts[i] = mu_left[i] - mu_right[i];

                // compute the estimate of mu for this site and node
                mu_node[i] = (mu_left[i] * t_right + mu_right[i] * t_left) / (t_left + t_right);
            }
            
            double lnl_contrast = RbStatistics::MultivariateNormal::lnPdfPrecision(means, precision_matrices[active_matrix], these_contrasts, branch_length);
            p_node = lnl_contrast + p_left + p_right;
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.

}


void PhyloMultivariateBrownianProcessREML::recursiveComputeContrasts( const TopologyNode &node, size_t node_index )
{
    
    // nothing to do for tips
    if ( node.isTip() == false )
    {

        std::vector<double> &mu_node = this->independent_contrasts[node_index];
        double              &sd_node = this->independent_contrasts_sds[node_index];

        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        
        for (size_t j = 1; j < num_children; ++j)
        {
            
            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left = &node.getChild(0);
                left_index = left->getIndex();
                recursiveComputeContrasts( *left, left_index );
            }
            
            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeContrasts( right, right_index );
            
            // get the per node and site contrasts
            const NodeCache &left_cache = this->node_likelihoods[left_index];
            const NodeCache &right_cache = this->node_likelihoods[right_index];

            const std::vector<double> &mu_left  = left_cache.contrasts;
            const std::vector<double> &mu_right = right_cache.contrasts;
            
            // get the propagated uncertainties
            double delta_left  = left_cache.contrast_uncertainty;
            double delta_right = right_cache.contrast_uncertainty;
            
            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());
            
            // add the propagated uncertainty to the branch lengths
            double t_left        = v_left  + delta_left;
            double t_right       = v_right + delta_right;
            double branch_length = t_left + t_right;
            
            sd_node = pow(branch_length, 0.5);
            for (size_t i = 0; i < this->num_sites; ++i)
            {
//                mu_node[i] = (mu_left[i] - mu_right[i]) / pow(branch_length, 0.5);
//                mu_node[i] = (mu_left[i] - mu_right[i]) / branch_length;
                mu_node[i] = mu_left[i] - mu_right[i];
            }
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}


void PhyloMultivariateBrownianProcessREML::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already invalid, then all ancestral nodes must have been invalidated too
    if ( node_likelihoods.is_valid(index) )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }

        node_likelihoods.invalidate(index);
        
    }
    
}


/*
 * Invalidate the computation that uses a branch transform.
 * Branch parameters belong to the edge ending at n, and are applied by n's parent.
 */
void PhyloMultivariateBrownianProcessREML::invalidateBranchAndAncestors( const TopologyNode &n )
{
    if (n.isRoot())
    {
        recursivelyFlagNodeDirty(n);
    }
    else
    {
        recursivelyFlagNodeDirty(n.getParent());
    }
}


/*
 * Invalidate recomputed multivariate Brownian REML entries while keeping fixed tip observations valid.
 * This is the full-cache fallback for topology and global parameter changes.
 */
void PhyloMultivariateBrownianProcessREML::invalidateInternalNodes( void )
{
    const std::vector<TopologyNode*> &nodes = this->tau->getValue().getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ((*it)->isTip() == false)
        {
            node_likelihoods.invalidate((*it)->getIndex());
        }
    }
}


void PhyloMultivariateBrownianProcessREML::resetValue( void )
{
    const bool had_snapshot = node_likelihoods.has_snapshot();

    this->num_nodes = tau->getValue().getNumberOfNodes();
    node_likelihoods.resize(this->num_nodes);
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = std::vector<size_t>(this->num_sites,0);
    size_t site_index = 0;
    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(site_index) )
        {
            ++site_index;
            if ( site_index >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }
        site_indices[i] = site_index;
        ++site_index;
    }
    
    const std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            size_t index = (*it)->getIndex();
            ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );

            NodeCache &tip_cache = node_likelihoods.init_for_writing(index);
            tip_cache.partial_likelihood = 0.0;
            tip_cache.contrasts.resize(this->num_sites);
            tip_cache.contrast_uncertainty = 0.0;

            for (size_t site = 0; site < this->num_sites; ++site)
            {
                double &c = taxon.getCharacter(site_indices[site]);
                tip_cache.contrasts[site] = c;
            }
        }
    }

    if (had_snapshot)
    {
        // Compatibility note: resetValue() can run during listener reattachment inside a proposal.
        // Replace resize()'s all-invalid rollback state with one that can recompute from valid tips.
        node_likelihoods.keep();
        node_likelihoods.snapshot();
    }

}


void PhyloMultivariateBrownianProcessREML::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the precision matrix if necessary
    if ( affecter == rate_matrix )
    {
        // Legacy local two-slot rollback for precision matrices.
        // Remove this once the matrix cache uses explicit snapshot/restore state.
        active_matrix = (active_matrix == 0 ? 1 : 0);
    }

    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.restore();
    }
}


void PhyloMultivariateBrownianProcessREML::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    
    // get the sequence of this node
    size_t node_index = node.getIndex();
    const ContinuousTaxonData &parent = taxa[ node_index ];
    
    std::vector<double> parent_state(num_sites, 0.0);
    for (size_t i = 0; i < num_sites; ++i)
    {
        parent_state[i] = parent.getCharacter(i);
    }
    
    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        
        // get the branch length for this child
        size_t child_index   = child.getIndex();
        double branch_length = this->computeBranchTime(child_index, child.getBranchLength());
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
//        std::vector<double> c = RbStatistics::MultivariateNormal::rvPrecision(parent_state, precision_matrices[active_matrix], *rng, branch_length);
        std::vector<double> c = RbStatistics::MultivariateNormal::rvCovariance(parent_state, rate_matrix->getValue(), *rng, branch_length);
        
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // add the character to the sequence
            taxon.addCharacter( c[i] );
        }
        
        if ( child.isTip() )
        {
            taxon.setTaxon( child.getTaxon() );
        }
        else
        {
            // recursively simulate the sequences
            simulateRecursively( child, taxa );
        }
        
    }
    
}


std::vector<double> PhyloMultivariateBrownianProcessREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


/*
 * Snapshot the per-node likelihood cache before invalidation mutates active cache metadata.
 * The derived independent-contrast vectors are recomputed on demand and are not snapshotted.
 */
void PhyloMultivariateBrownianProcessREML::snapshotSpecialization( void )
{
    node_likelihoods.snapshot();
}


/*
 * Mark multivariate Brownian REML likelihood caches dirty after a dependency changes.
 * Rate-matrix invalidation also refreshes the active precision matrix.
 */
void PhyloMultivariateBrownianProcessREML::invalidateSpecialization( const DagNode* affecter, bool touchAll )
{
 
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_clock_rates )
    {
        
        const std::set<size_t> &indices = this->heterogeneous_clock_rates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everything for recomputation
            touchAll = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                if ( *it >= nodes.size() )
                {
                    touchAll = true;
                    break;
                }

                this->invalidateBranchAndAncestors( *nodes[*it] );
            }
        }
    }
    else if ( affecter == rate_matrix )
    {
        // compute the inverse variance-covariance matrix (the precision matrix)
        active_matrix = (active_matrix == 0 ? 1 : 0);
        precision_matrices[active_matrix] = rate_matrix->getValue().computeInverse();
        precision_matrices[active_matrix].setCholesky(true);
        
        // we need to recompute the likelihood
        touchAll = true;
    }
    else if ( affecter == static_cast<const DagNode*>(this->dag_node) )
    {
        resetValue();
    }
    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
    }
    
    if ( touchAll )
    {
        invalidateInternalNodes();
    }
    
}


/** Swap a parameter of the distribution */
void PhyloMultivariateBrownianProcessREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    if (oldP == this->rate_matrix)
    {
        rate_matrix = static_cast<const TypedDagNode< MatrixReal >* >( newP );
        active_matrix = 0;
        precision_matrices[0] = rate_matrix->getValue().computeInverse();
        precision_matrices[0].setCholesky(true);
        precision_matrices[1] = rate_matrix->getValue().computeInverse();
        precision_matrices[1].setCholesky(true);
    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}
