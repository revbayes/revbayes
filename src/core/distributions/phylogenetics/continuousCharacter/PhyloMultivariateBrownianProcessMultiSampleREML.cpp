#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "DistributionNormal.h"
#include "DistributionMultivariateNormal.h"
#include "PhyloMultivariateBrownianProcessMultiSampleREML.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloMultivariateBrownianProcessMultiSampleREML::PhyloMultivariateBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode<MatrixReal> *c, const TypedDagNode< RbVector< RbVector< double > > > *v, const std::vector<Taxon> &ta, size_t ns) : AbstractPhyloBrownianProcess( tr, ns ),
//PhyloMultivariateBrownianProcessMultiSampleREML::PhyloMultivariateBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode<MatrixReal> *c, const TypedDagNode< RbVector< double > > *v, const TypedDagNode< RbVector< double > > *v2, const std::vector<Taxon> &ta, size_t ns) : AbstractPhyloBrownianProcess( tr, ns ),
    within_species_variances( v ),
//    within_species_variances2( v2 ),
    node_likelihoods(this->num_nodes),
    taxa( ta ),
    rate_matrix( c )
{
    
    // add the parameters to our set
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( within_species_variances );
//    this->addParameter( within_species_variances2 );
    this->addParameter( rate_matrix );
    
    // make sure the rate matrix is inverted using Cholesky decomposition
    rate_matrix->getValue().setCholesky(true);
    
    // We don'e want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // now we need to reset the value
    this->redrawValue();
    
    // we need to reset the contrasts
    resetValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloMultivariateBrownianProcessMultiSampleREML::~PhyloMultivariateBrownianProcessMultiSampleREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloMultivariateBrownianProcessMultiSampleREML* PhyloMultivariateBrownianProcessMultiSampleREML::clone( void ) const
{
    
    return new PhyloMultivariateBrownianProcessMultiSampleREML( *this );
}


double PhyloMultivariateBrownianProcessMultiSampleREML::computeLnProbability( void )
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
    size_t root_index = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( not this->node_likelihoods.is_valid(root_index) )
    {
        recursiveComputeLnProbability( root, root_index );
    }

    // sum the partials up
    this->ln_prob = this->node_likelihoods[root_index].partial_likelihood;

    return this->ln_prob;
}


double PhyloMultivariateBrownianProcessMultiSampleREML::computeMeanForSpecies(const std::string &name, size_t index)
{
    
    double mean = 0.0;
    double num_samples = 0.0;
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        
        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
            
            if ( taxon.isCharacterResolved( index ) == true )
            {
                mean += taxon.getCharacter(index);
            
                ++num_samples;
            }
            
        }
        
    }
    
    // normalize
    mean /= num_samples;
    
    
    return mean;
}


void PhyloMultivariateBrownianProcessMultiSampleREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


double PhyloMultivariateBrownianProcessMultiSampleREML::getNumberOfSamplesForSpecies(const std::string &name)
{
    
    double num_samples = 0.0;
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        
        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ++num_samples;
        }
        
    }
    
    return num_samples;
}


double PhyloMultivariateBrownianProcessMultiSampleREML::getWithinSpeciesVariance(const std::string &name, size_t site)
{
    
    size_t index = this->tau->getValue().getTipIndex( name );

//    if (site == 0)
//    {
//        return within_species_variances->getValue()[ index ];
//    }
//    else
//    {
//        return within_species_variances2->getValue()[ index ];
//    }
    return within_species_variances->getValue()[site][ index ];
}




void PhyloMultivariateBrownianProcessMultiSampleREML::keepSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.keep();
    }
}


void PhyloMultivariateBrownianProcessMultiSampleREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    // check for recomputation
    if ( node.isTip() == true && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihood = 0.0;
        node_cache.contrasts.assign(this->num_sites, 0.0);
        node_cache.contrast_uncertainty.assign(this->num_sites, 0.0);

        double &p_node  = node_cache.partial_likelihood;
        std::vector<double> &mu_node  = node_cache.contrasts;
        
        const std::string &name = this->tau->getValue().getNode( node_index ).getName();
        double num_samples = 0.0;

        for (int i=0; i<this->num_sites; i++)
        {
            mu_node[i] = computeMeanForSpecies(name, site_indices[i]);
        }

        for (size_t i=0; i<taxa.size(); ++i)
        {
            
            const Taxon &t = taxa[i];
            if ( name == t.getSpeciesName() )
            {
                
                for (int i=0; i<this->num_sites; i++)
                {
                    double var = getWithinSpeciesVariance(name, i);
                    double stdev = sqrt( var );

                    ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
                    
                    if ( taxon.isCharacterResolved( site_indices[i] ) == true )
                    {
                        double x = taxon.getCharacter( site_indices[i] );
                    
                        // get the site specific rate of evolution
                        double standDev = stdev;
                    
                        // compute the contrasts for this site and node
                        double contrast = mu_node[i] - x;
                    
                        // compute the probability for the contrasts at this node
                        double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);
                    
                        // sum up the probabilities of the contrasts
                        p_node += lnl_node;
                    
                    } // if is resolved
                    
                } // end for-loop over all sites
                
                
                ++num_samples;
            }
            
        }
        
        for (int site=0; site<this->num_sites; site++)
        {
            double var = getWithinSpeciesVariance(name, site);
            node_cache.contrast_uncertainty[site] = var / num_samples;
        }
        
    }
    else if ( node.isTip() == false && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihood = 0.0;
        node_cache.contrasts.assign(this->num_sites, 0.0);
        node_cache.contrast_uncertainty.assign(this->num_sites, 0.0);

        double &p_node  = node_cache.partial_likelihood;
        std::vector<double> &mu_node  = node_cache.contrasts;
        
        
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
            
            double p_left  = left_cache->partial_likelihood;
            double p_right = right_cache.partial_likelihood;
            
            // get the per node and site contrasts
            const std::vector<double> &mu_left  = left_cache->contrasts;
            const std::vector<double> &mu_right = right_cache.contrasts;
            
            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());
            
            std::vector<double> these_contrasts(num_sites, 0.0);
            std::vector<double> means(num_sites, 0.0);
            std::vector<double> these_branch_lengths(num_sites, 0.0);
            
//            double stdev = sqrt(t_left+t_right);
            for (int i=0; i<this->num_sites; i++)
            {
                // get the propagated uncertainties
                double delta_left  = left_cache->contrast_uncertainty[i];
                double delta_right = right_cache.contrast_uncertainty[i];
                
                // add the propagated uncertainty to the branch lengths
                double t_left  = v_left  + delta_left;
                double t_right = v_right + delta_right;
                
                these_branch_lengths[i] = sqrt(t_left + t_right);

                // set delta_node = (t_l*t_r)/(t_l+t_r);
                node_cache.contrast_uncertainty[i] = (t_left*t_right) / (t_left+t_right);

                // compute the contrasts for this site and node
                these_contrasts[i] = mu_left[i] - mu_right[i];
                
                // compute the estimate of mu for this site and node
                mu_node[i] = (mu_left[i]*t_right + mu_right[i]*t_left) / (t_left+t_right);

            } // end for-loop over all sites
            
            double lnl_contrast = RbStatistics::MultivariateNormal::lnPdfCovariance(means, rate_matrix->getValue(), these_contrasts, these_branch_lengths);
//            double lnl_contrast = RbStatistics::MultivariateNormal::lnPdfCovariance(means, rate_matrix->getValue(), these_contrasts, 1.0);
            p_node = lnl_contrast + p_left + p_right;

        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}


void PhyloMultivariateBrownianProcessMultiSampleREML::recursivelyFlagNodeDirty( const TopologyNode &n )
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
void PhyloMultivariateBrownianProcessMultiSampleREML::invalidateBranchAndAncestors( const TopologyNode &n )
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



void PhyloMultivariateBrownianProcessMultiSampleREML::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new ContinuousCharacterData();
    
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxon_data = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    // simulate the root sequence
    ContinuousTaxonData &root = taxon_data[ tau->getValue().getRoot().getIndex() ];
    
    std::vector<double> root_states = simulateRootCharacters(num_sites);
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // create the character
        double c = root_states[i];
        
        // add the character to the sequence
        root.addCharacter( c );
    }
    
    // recursively simulate the sequences
    simulateRecursively( tau->getValue().getRoot(), taxon_data );
    
    // Get the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        const std::string &species_name = tau->getValue().getNode(i).getName();
        const ContinuousTaxonData &species_data = taxon_data[i];
        
        for ( size_t j=0; j<taxa.size(); ++j )
        {
            
            const Taxon &t = taxa[j];
            if ( species_name == t.getSpeciesName() )
            {
                ContinuousTaxonData individual_data = ContinuousTaxonData( t );
                
                for ( size_t k = 0; k < num_sites; ++k )
                {
                    double species_sigma = sqrt( getWithinSpeciesVariance( species_name, k ) );
                    
                    // get the ancestral character for this site
                    double parent_state = species_data.getCharacter(k);
                    
                    // compute the standard deviation for this site
                    double stand_dev = species_sigma * computeSiteRate( k );
                    
                    // create the character
                    double c = RbStatistics::Normal::rv( parent_state, stand_dev, *rng);
                    
                    // add the character to the sequence
                    individual_data.addCharacter( c );
                }
                
                this->value->addTaxonData( individual_data );
            }
            
        }
    }
    
    // tell the derived classes
    this->resetValue();
}


void PhyloMultivariateBrownianProcessMultiSampleREML::resetValue( void )
{
    this->num_nodes = tau->getValue().getNumberOfNodes();
    node_likelihoods.resize(this->num_nodes);

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    site_indices = std::vector<size_t>(this->num_sites,0);
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
    
}


void PhyloMultivariateBrownianProcessMultiSampleREML::restoreSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.restore();
    }
}


std::vector<double> PhyloMultivariateBrownianProcessMultiSampleREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


void PhyloMultivariateBrownianProcessMultiSampleREML::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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


/*
 * Snapshot the per-node likelihood cache before invalidation mutates active cache metadata.
 * Tips are recomputed cache entries because within-species variance affects their likelihood.
 */
void PhyloMultivariateBrownianProcessMultiSampleREML::snapshotSpecialization( void )
{
    node_likelihoods.snapshot();
}


/*
 * Mark multisample multivariate Brownian REML caches dirty after a dependency changes.
 * Rate-matrix invalidation marks everything dirty because likelihoods read the covariance directly.
 */
void PhyloMultivariateBrownianProcessMultiSampleREML::invalidateSpecialization( const DagNode* affecter, bool touchAll )
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
    else if ( affecter == rate_matrix )
    {
        // we need to recompute the likelihood
        touchAll = true;
    }
    else if ( affecter == this->within_species_variances )
    {
        // Compatibility note: nested vector touched indices do not identify both site and species.
        // Keep full invalidation until nested touched-index metadata can address individual tips.
        touchAll = true;
        
    }
//    else if ( affecter == this->within_species_variances2 )
//    {
//
//        const std::set<size_t> &indices = this->within_species_variances2->getTouchedElementIndices();
//
//        // maybe all of them have been touched or the flags haven't been set properly
//        if ( indices.size() == 0 )
//        {
//            // just flag everyting for recomputation
//            touchAll = true;
//        }
//        else
//        {
//            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
//            // flag recomputation only for the nodes
//            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
//            {
//                this->recursivelyFlagNodeDirty( *nodes[*it] );
//
//            }
//        }
//
//    }
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
        node_likelihoods.invalidate_all();
    }
    
}


/** Swap a parameter of the distribution */
void PhyloMultivariateBrownianProcessMultiSampleREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else if (oldP == this->rate_matrix)
    {
        rate_matrix = static_cast<const TypedDagNode< MatrixReal >* >( newP );
    }
    else if ( oldP == this->within_species_variances )
    {
//        within_species_variances = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
        within_species_variances = static_cast<const TypedDagNode< RbVector< RbVector< double > > >* >( newP );
    }
//    else if ( oldP == this->within_species_variances2 )
//    {
//        within_species_variances2 = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
//    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}
