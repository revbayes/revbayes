#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "DistributionNormal.h"
#include "PhyloBrownianProcessREML.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

PhyloBrownianProcessREML::PhyloBrownianProcessREML(const TypedDagNode<Tree> *t, size_t ns) :
    AbstractPhyloBrownianProcess( t, ns ),
    node_likelihoods(this->num_nodes)
{
    
    
    // We don'e want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // now we need to reset the value
    this->redrawValue();
    
    // we need to reset the means
    resetValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloBrownianProcessREML::~PhyloBrownianProcessREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloBrownianProcessREML* PhyloBrownianProcessREML::clone( void ) const
{
    
    return new PhyloBrownianProcessREML( *this );
}


double PhyloBrownianProcessREML::computeLnProbability( void )
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
        
        // start by filling the likelihood vector for the children of the root
        if ( root.getNumberOfChildren() != 2 && root.getNumberOfChildren() != 3 ) // rooted trees have two children for the root
        {
            throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
        }
    }

    // NOTE: After restoreSpecialization( ) is called, all the nodes will be marked clean.
    //       But we still need to update ln_prob

    // sum the partials up
    ln_prob = sumRootLikelihood();

    return ln_prob;
}



void PhyloBrownianProcessREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


void PhyloBrownianProcessREML::keepSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.keep();
    }
}


void PhyloBrownianProcessREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{

    // check for recomputation
    if ( node.isTip() == false && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihoods.assign(this->num_sites, 0.0);
        node_cache.means.assign(this->num_sites, 0.0);
        node_cache.variances_per_site.assign(this->num_sites, 0.0);
        node_cache.missing_data.assign(this->num_sites, false);
        node_cache.variance = 0.0;

        std::vector<double> &p_node   = node_cache.partial_likelihoods;
        std::vector<double> &mu_node  = node_cache.means;

        
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
            const std::vector<double> &p_left  = left_cache->partial_likelihoods;
            const std::vector<double> &p_right = right_cache.partial_likelihoods;

            // get the per node and site means
            const std::vector<double> &mu_left  = left_cache->means;
            const std::vector<double> &mu_right = right_cache.means;
            
            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());
            
            // get the propagated uncertainties
            double delta_left   = 0.0;
            double delta_right  = 0.0;
            double var_left     = 0.0;
            double var_right    = 0.0;
            double stdev        = 0.0;
            if ( use_missing_data == false )
            {
                delta_left  = left_cache->variance;
                delta_right = right_cache.variance;

                // add the propagated uncertainty to the branch lengths
                var_left  = v_left  + delta_left;
                var_right = v_right + delta_right;

                // set delta_node = (t_l*t_r)/(t_l+t_r);
                node_cache.variance = (var_left*var_right) / (var_left+var_right);

                stdev = sqrt(var_left+var_right);
            }

            for (int i=0; i<this->num_sites; ++i)
            {

                if ( use_missing_data == true )
                {
                    delta_left  = left_cache->variances_per_site[i];
                    delta_right = right_cache.variances_per_site[i];

                    // add the propagated uncertainty to the branch lengths
                    var_left  = v_left  + delta_left;
                    var_right = v_right + delta_right;

                    // set delta_node = (t_l*t_r)/(t_l+t_r);
                    node_cache.variances_per_site[i] = (var_left*var_right) / (var_left+var_right);
                    
                    stdev = sqrt(var_left+var_right);
                }
                
                if ( use_missing_data == true && left_cache->missing_data[i] == true && right_cache.missing_data[i] == true )
                {
                    node_cache.missing_data[i] = true;
                    
                    p_node[i]  = p_left[i] + p_right[i];
                    mu_node[i] = RbConstants::Double::nan;

                    node_cache.variances_per_site[i] = 0.0;
                }
                else if ( use_missing_data == true && left_cache->missing_data[i] == true && right_cache.missing_data[i] == false )
                {
                    node_cache.missing_data[i] = false;
                    
                    p_node[i]  = p_left[i] + p_right[i];
                    mu_node[i] = mu_right[i];
                    
                    node_cache.variances_per_site[i] = var_right;

                }
                else if ( use_missing_data == true && left_cache->missing_data[i] == false && right_cache.missing_data[i] == true )
                {
                    node_cache.missing_data[i] = false;
                    
                    p_node[i]  = p_left[i] + p_right[i];
                    mu_node[i] = mu_left[i];
                    
                    node_cache.variances_per_site[i] = var_left;
                }
                else
                {

                    // get the site specific rate of evolution
                    double standDev = this->computeSiteRate(i) * stdev;

                    // compute the means for this site and node
                    double contrast = mu_left[i] - mu_right[i];

                    // compute the probability for the means at this node
                    double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);

                    // sum up the probabilities of the means
                    p_node[i] = lnl_node + p_left[i] + p_right[i];
                    
                    mu_node[i] = (mu_left[i]*var_right + mu_right[i]*var_left) / (var_left+var_right);
                    
                    if ( use_missing_data == true )
                    {
                        node_cache.missing_data[i] = false;
                        node_cache.variances_per_site[i] = (var_left*var_right) / (var_left+var_right);
                    }
                    
                }
                

            } // end for-loop over all sites

        } // end for-loop over all children
        
    } // end if we need to compute something for this node.

}



void PhyloBrownianProcessREML::recursivelyFlagNodeDirty( const TopologyNode &n )
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
void PhyloBrownianProcessREML::invalidateBranchAndAncestors( const TopologyNode &n )
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
 * Invalidate all recomputed Brownian REML entries while leaving fixed tip observations clean.
 * This is the full-cache fallback for topology and global parameter changes.
 */
void PhyloBrownianProcessREML::invalidateInternalNodes( void )
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


void PhyloBrownianProcessREML::resetValue( void )
{
    
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
    
    // first we check for missing data
    use_missing_data = false;
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();

    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            size_t index = (*it)->getIndex();
            ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );

            NodeCache &tip_cache = node_likelihoods.init_for_writing(index);
            tip_cache.partial_likelihoods.assign(this->num_sites, 0.0);
            tip_cache.means.resize(this->num_sites);
            tip_cache.variances_per_site.assign(this->num_sites, 0.0);
            tip_cache.missing_data.assign(this->num_sites, false);
            tip_cache.variance = 0.0;

            for (size_t site = 0; site < this->num_sites; ++site)
            {
                double c = taxon.getCharacter(site_indices[site]);

                tip_cache.means[site] = c;

                if ( RbMath::isFinite(c) == false )
                {
                    tip_cache.missing_data[site] = true;
                    use_missing_data = true;
                }
            }
        }
    }
    
}


void PhyloBrownianProcessREML::restoreSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.restore();
    }
}


std::vector<double> PhyloBrownianProcessREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


double PhyloBrownianProcessREML::sumRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // get the index of the root node
    size_t node_index = root.getIndex();
    
    // get the pointers to the partial likelihoods of the left and right subtree
    const std::vector<double> &p_node = this->node_likelihoods[node_index].partial_likelihoods;
    
    // sum the log-likelihoods for all sites together
    double sum_partial_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        sum_partial_probs += p_node[site];
    }
    
    return sum_partial_probs;
}


void PhyloBrownianProcessREML::snapshotSpecialization( void )
{
    node_likelihoods.snapshot();
}


/*
 * Mark Brownian REML likelihood caches dirty after a dependency changes.
 * Snapshot state is stored by IndexedSnapshotCache before this invalidation hook runs.
 */
void PhyloBrownianProcessREML::invalidateSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_clock_rates )
    {
        
        const std::set<size_t> &indices = this->heterogeneous_clock_rates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everyting for recomputation
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

    if ( affecter == this->dag_node )
    {
        resetValue();
    }
    else if ( affecter != this->heterogeneous_clock_rates && affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
    }
    
    if ( touchAll )
    {
        invalidateInternalNodes();
    }
    
}


/** Swap a parameter of the distribution */
void PhyloBrownianProcessREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}
