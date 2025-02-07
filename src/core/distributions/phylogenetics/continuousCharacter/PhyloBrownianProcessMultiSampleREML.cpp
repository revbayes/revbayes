#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "DistributionNormal.h"
#include "PhyloBrownianProcessMultiSampleREML.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloBrownianProcessMultiSampleREML::PhyloBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode< RbVector< double > > *v, const std::vector<Taxon> &ta, size_t ns) : AbstractPhyloBrownianProcess( tr, ns ),
    within_species_variances( v ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrasts( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrast_uncertainty( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
    taxa( ta )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( within_species_variances );
    
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
PhyloBrownianProcessMultiSampleREML::~PhyloBrownianProcessMultiSampleREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloBrownianProcessMultiSampleREML* PhyloBrownianProcessMultiSampleREML::clone( void ) const
{
    
    return new PhyloBrownianProcessMultiSampleREML( *this );
}


double PhyloBrownianProcessMultiSampleREML::computeLnProbability( void )
{
    
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        dirty_nodes = std::vector<bool>(tau->getValue().getNumberOfNodes(), true);
    }
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( this->dirty_nodes[root_index] )
        recursiveComputeLnProbability( root, root_index );

    // sum the partials up
    this->ln_prob = sumRootLikelihood();

    return this->ln_prob;
}


double PhyloBrownianProcessMultiSampleREML::computeMeanForSpecies(const std::string &name, size_t index)
{
    
    double mean = 0.0;
    double num_samples = 0.0;
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        
        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
            
            if ( taxon.isCharacterResolved(index) == true )
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


void PhyloBrownianProcessMultiSampleREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


double PhyloBrownianProcessMultiSampleREML::getNumberOfSamplesForSpecies(const std::string &name)
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


double PhyloBrownianProcessMultiSampleREML::getWithinSpeciesVariance(const std::string &name)
{
    
    size_t index = this->tau->getValue().getTipIndex( name );
    
    return within_species_variances->getValue()[ index ];
}




void PhyloBrownianProcessMultiSampleREML::keepSpecialization( const DagNode* affecter )
{
    
    // reset all flags
    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    for (std::vector<bool>::iterator it = this->changed_nodes.begin(); it != this->changed_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
}


void PhyloBrownianProcessMultiSampleREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    // check for recomputation
    if ( node.isTip() == true && dirty_nodes[node_index] == true )
    {
        
        std::vector<double> &p_node  = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node  = this->contrasts[this->active_likelihood[node_index]][node_index];

        const std::string &name = this->tau->getValue().getNode( node_index ).getName();
        double num_samples = 0.0;
        
        double var = getWithinSpeciesVariance(name);
        
        double stdev = sqrt( var );
        
        for (int i=0; i<this->num_sites; i++)
        {
             p_node[i] = 0;
        }
        
        for (size_t i=0; i<taxa.size(); ++i)
        {
            
            const Taxon &t = taxa[i];
            if ( name == t.getSpeciesName() )
            {
                
                for (int i=0; i<this->num_sites; i++)
                {
                    
                    ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
                    if ( taxon.isCharacterResolved( site_indices[i] ) )
                    {
                        double x = taxon.getCharacter( site_indices[i] );
                    
                        // get the site specific rate of evolution
                        double standDev = this->computeSiteRate(i) * stdev;
                    
                        // compute the contrasts for this site and node
                        double contrast = mu_node[i] - x;
                    
                        // compute the probability for the contrasts at this node
                        double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);
                    
                        // sum up the probabilities of the contrasts
                        p_node[i] += lnl_node;
                        
                    } // end if is character resolved
                    
                } // end for-loop over all sites

                
                ++num_samples;
            }
            
        }
        
        contrast_uncertainty[this->active_likelihood[node_index]][node_index] = var / num_samples;
        
    }
    else if ( node.isTip() == false && dirty_nodes[node_index] == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        std::vector<double> &p_node  = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node  = this->contrasts[this->active_likelihood[node_index]][node_index];
        
        
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
            
            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeLnProbability( right, right_index );
            
            const std::vector<double> &p_left  = this->partial_likelihoods[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &p_right = this->partial_likelihoods[this->active_likelihood[right_index]][right_index];
            
            // get the per node and site contrasts
            const std::vector<double> &mu_left  = this->contrasts[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->contrasts[this->active_likelihood[right_index]][right_index];
            
            // get the propagated uncertainties
            double delta_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
            double delta_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];
            
            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, tau->getValue().getBranchLengthForNode(*left));
            }
            double v_right = this->computeBranchTime(right_index, tau->getValue().getBranchLengthForNode(right));
            
            // add the propagated uncertainty to the branch lengths
            double t_left  = v_left  + delta_left;
            double t_right = v_right + delta_right;
            
            // set delta_node = (t_l*t_r)/(t_l+t_r);
            this->contrast_uncertainty[this->active_likelihood[node_index]][node_index] = (t_left*t_right) / (t_left+t_right);
            
            double stdev = sqrt(t_left+t_right);
            for (int i=0; i<this->num_sites; i++)
            {
                
                // get the site specific rate of evolution
                double standDev = this->computeSiteRate(i) * stdev;
                
                // compute the contrasts for this site and node
                double contrast = mu_left[i] - mu_right[i];
                
                // compute the probability for the contrasts at this node
                double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);
                
                // sum up the probabilities of the contrasts
                p_node[i] = lnl_node + p_left[i] + p_right[i];
                
                // compute the estimate of mu for this site and node
                mu_node[i] = (mu_left[i]*t_right + mu_right[i]*t_left) / (t_left+t_right);

            } // end for-loop over all sites
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}


void PhyloBrownianProcessMultiSampleREML::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( dirty_nodes[index] == false )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }
        
        // set the flag
        dirty_nodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_nodes[index] == false )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }
        
    }
    
}



void PhyloBrownianProcessMultiSampleREML::redrawValue( void )
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
        double species_sigma = sqrt( getWithinSpeciesVariance( species_name ) );
        
        for ( size_t j=0; j<taxa.size(); ++j )
        {
            
            const Taxon &t = taxa[j];
            if ( species_name == t.getSpeciesName() )
            {
                ContinuousTaxonData individual_data = ContinuousTaxonData( t );
                
                for ( size_t k = 0; k < num_sites; ++k )
                {
                    
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


void PhyloBrownianProcessMultiSampleREML::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrasts = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrast_uncertainty = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );
    
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
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                const std::string &name = (*it)->getName();
                double c = computeMeanForSpecies(name, site_indices[site]);
                contrasts[0][(*it)->getIndex()][site] = c;
                contrasts[1][(*it)->getIndex()][site] = c;
                contrast_uncertainty[0][(*it)->getIndex()] = sqrt( getWithinSpeciesVariance(name) ) / getNumberOfSamplesForSpecies(name);
                contrast_uncertainty[1][(*it)->getIndex()] = sqrt( getWithinSpeciesVariance(name) ) / getNumberOfSamplesForSpecies(name);
            }
        }
    }
    
    
    // finally we set all the flags for recomputation
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = true;
    }
    
    // flip the active likelihood pointers
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        active_likelihood[index] = 0;
        changed_nodes[index] = true;
    }
    
}


void PhyloBrownianProcessMultiSampleREML::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    // restore the active likelihoods vector
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_nodes[index] == true )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
        }
        
        // set all flags to false
        changed_nodes[index] = false;
    }
    
}


std::vector<double> PhyloBrownianProcessMultiSampleREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


double PhyloBrownianProcessMultiSampleREML::sumRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // get the index of the root node
    size_t node_index = root.getIndex();
    
    // get the pointers to the partial likelihoods of the left and right subtree
    std::vector<double> &p_node = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
    
    // sum the log-likelihoods for all sites together
    double sum_partial_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        sum_partial_probs += p_node[site];
    }
    
    return sum_partial_probs;
}


void PhyloBrownianProcessMultiSampleREML::touchSpecialization( const DagNode* affecter, bool touchAll )
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
                if ( nodes[*it]->isRoot() == false )
                {
                    this->recursivelyFlagNodeDirty( nodes[*it]->getParent() );
                }
            }
        }
    }
    else if ( affecter == this->within_species_variances )
    {
        
        const std::set<size_t> &indices = this->within_species_variances->getTouchedElementIndices();
        
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
                this->recursivelyFlagNodeDirty( *nodes[*it] );
                
            }
        }
        
    }
    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
        
        if ( affecter == this->dag_node )
        {
            resetValue();
        }
        
    }
    
    if ( touchAll )
    {
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }
        
        // flip the active likelihood pointers
        for (size_t index = 0; index < changed_nodes.size(); ++index)
        {
            if ( changed_nodes[index] == false )
            {
                active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
                changed_nodes[index] = true;
            }
        }
        
    }
    
}


/** Swap a parameter of the distribution */
void PhyloBrownianProcessMultiSampleREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else if ( oldP == this->within_species_variances )
    {
        within_species_variances = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}

