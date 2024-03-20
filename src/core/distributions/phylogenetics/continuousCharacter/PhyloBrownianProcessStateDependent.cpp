#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "AbstractPhyloContinuousCharacterProcess.h"
#include "BranchHistory.h"
#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "PhyloBrownianProcessStateDependent.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StandardState.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TreeHistoryCtmc.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloBrownianProcessStateDependent::PhyloBrownianProcessStateDependent(const TypedDagNode<CharacterHistoryDiscrete> *ch, size_t ns) : TypedDistribution< ContinuousCharacterData > ( new ContinuousCharacterData() ),
    num_nodes( ch->getValue().getNumberBranches()+1 ),
    num_sites( ns ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrasts( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrast_uncertainty( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    contrast_uncertainty_per_site( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
    character_histories( ch ),
    use_missing_data(false)
{
    // initialize default parameters
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    state_dependent_sigma       = NULL;
    
    
    // add parameters
    addParameter( homogeneous_sigma );
    addParameter( character_histories );
    
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
PhyloBrownianProcessStateDependent::~PhyloBrownianProcessStateDependent( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
//    // remove myself from the tree listeners
//    if ( tau != NULL )
//    {
//        tau->getValue().getTreeChangeEventHandler().removeListener( this );
//    }
    
}



PhyloBrownianProcessStateDependent* PhyloBrownianProcessStateDependent::clone( void ) const
{
    
    return new PhyloBrownianProcessStateDependent( *this );
}



double PhyloBrownianProcessStateDependent::computeBranchTime( size_t node_idx, double brlen )
{
    
    // get the clock rate for the branch
    double branch_time = 1.0;
//    if ( this->heterogeneous_clock_rates != NULL )
//    {
//        double sigma = this->heterogeneous_clock_rates->getValue()[node_idx];
//        branch_time = sigma * sigma * brlen;
//    }
//    else
//    {
//        double sigma = this->homogeneous_clock_rate->getValue();
//        branch_time = sigma * sigma * brlen;
//    }
    
    return branch_time;
}


double PhyloBrownianProcessStateDependent::computeSiteRate( size_t site_idx )
{
    
    double rate = 1.0;
    
    return rate;
}




double PhyloBrownianProcessStateDependent::computeStateDependentSigma(size_t state_idx) const
{
    
    // get the drift rate for the branch
    double s;
    if ( this->state_dependent_sigma != NULL )
    {
        s = this->state_dependent_sigma->getValue()[state_idx];
    }
    else
    {
        s = this->homogeneous_sigma->getValue();
    }
    
    return s;
}


double PhyloBrownianProcessStateDependent::computeLnProbability( void )
{
    
//    // we need to check here if we still are listining to this tree for change events
//    // the tree could have been replaced without telling us
//    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
//    {
//        tau->getValue().getTreeChangeEventHandler().addListener( this );
//        dirty_nodes = std::vector<bool>(tau->getValue().getNumberOfNodes(), true);
//    }
    
    const Tree& tau = character_histories->getValue().getTree();
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = tau.getRoot();
    
    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
        
    // only necessary if the root is actually dirty
    if ( this->dirty_nodes[root_index] )
    {
                
        recursiveComputeLnProbability( root, root_index );
        
        // sum the partials up
        this->ln_prob = sumRootLikelihood();
        
    }
    return this->ln_prob;
}



//void PhyloBrownianProcessStateDependent::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
//{
//
//    // call a recursive flagging of all node above (closer to the root) and including this node
//    recursivelyFlagNodeDirty( n );
//
//}


void PhyloBrownianProcessStateDependent::keepSpecialization( const DagNode* affecter )
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


void PhyloBrownianProcessStateDependent::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    bool use_missing_data = false;
    
    // check for recomputation
    if ( node.isTip() == false && dirty_nodes[node_index] == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        std::vector<double> &p_node             = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node            = this->contrasts[this->active_likelihood[node_index]][node_index];
        
        
        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        
        for (size_t j = 1; j < num_children; ++j)
        {
            
            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left        = &node.getChild(0);
                left_index  = left->getIndex();
                recursiveComputeLnProbability( *left, left_index );
            }
                        
            const TopologyNode &right   = node.getChild(j);
            size_t right_index          = right.getIndex();
            recursiveComputeLnProbability( right, right_index );
            
            const std::vector<double> &p_left  = this->partial_likelihoods[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &p_right = this->partial_likelihoods[this->active_likelihood[right_index]][right_index];
            
            // get the per node and site contrasts
            const std::vector<double> &mu_left  = this->contrasts[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->contrasts[this->active_likelihood[right_index]][right_index];
            
            // get character history for tree
            const CharacterHistory& current_history = character_histories->getValue();

            // get the scaled branch lengths
            double v_left  = 0;
            // only compute the left branch if we didn't already
            if ( j == 1 )
            {
                double begin_time_left = left->getAge();
                
                // get branch history for node
                const BranchHistory& bh_left = current_history.getHistory(left_index);
                const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_left = bh_left.getHistory();
                
                for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_left.begin(); it!=hist_left.end(); ++it)
                {
                    CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                    double event_time = event->getAge();
                    
                    // we need to set the current rate category
                    size_t current_state = event->getState();
                    
                    double time_left = event_time - begin_time_left;
                    double sigma_left = computeStateDependentSigma( current_state );
                    
                    v_left += (sigma_left*sigma_left) * time_left;
                    
                    // update our variables
                    begin_time_left     = event_time;
                }
                size_t last_state_left  = static_cast<CharacterEventDiscrete*>(bh_left.getParentCharacters()[0])->getState();
                double sigma_left = computeStateDependentSigma( last_state_left );
                v_left += (sigma_left*sigma_left) * (node.getAge() - begin_time_left);
            }
            
            // get the scaled branch lengths (right side)
            double v_right  = 0;
            double begin_time_right = right.getAge();
            
            // get branch history for node
            const BranchHistory& bh_right = current_history.getHistory(right_index);
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_right = bh_right.getHistory();

            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_right.begin(); it!=hist_right.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                
                // we need to set the current rate category
                size_t current_state = event->getState();
                
                double time_right = event_time - begin_time_right;
                double sigma_right = computeStateDependentSigma( current_state );
                
                v_right += (sigma_right*sigma_right) * time_right;
                
                // update our variables
                begin_time_right     = event_time;
            }
            size_t last_state_right  = static_cast<CharacterEventDiscrete*>(bh_right.getParentCharacters()[0])->getState();
            double sigma_right = computeStateDependentSigma( last_state_right );
            v_right += (sigma_right*sigma_right) * (node.getAge() - begin_time_right);

            // get the propagated uncertainties
            double delta_left   = 0.0;
            double delta_right  = 0.0;
            double t_left       = 0.0;
            double t_right      = 0.0;
            double stdev        = 0.0;
            if ( use_missing_data == false )
            {
                delta_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
                delta_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];

                // add the propagated uncertainty to the branch lengths
                t_left  = v_left  + delta_left;
                t_right = v_right + delta_right;

                // set delta_node = (t_l*t_r)/(t_l+t_r);
                this->contrast_uncertainty[this->active_likelihood[node_index]][node_index] = (t_left*t_right) / (t_left+t_right);

                stdev = sqrt(t_left+t_right);
            }

            for (int site=0; site<this->num_sites; ++site)
            {

                if ( use_missing_data == true )
                {
//                    delta_left  = this->contrast_uncertainty_per_site[this->active_likelihood[left_index]][left_index][site];
//                    delta_right = this->contrast_uncertainty_per_site[this->active_likelihood[right_index]][right_index][site];

                    // add the propagated uncertainty to the branch lengths
                    t_left  = v_left  + delta_left;
                    t_right = v_right + delta_right;

                    // set delta_node = (t_l*t_r)/(t_l+t_r);
//                    this->contrast_uncertainty_per_site[this->active_likelihood[node_index]][node_index][site] = (t_left*t_right) / (t_left+t_right);
                    
                    stdev = sqrt(t_left+t_right);
                }
                
                if ( use_missing_data == true && missing_data[left_index][site] == true && missing_data[right_index][site] == true )
                {
                    missing_data[node_index][site] = true;
                    
                    p_node[site]  = p_left[site] + p_right[site];
                    mu_node[site] = RbConstants::Double::nan;

                    this->contrast_uncertainty_per_site[this->active_likelihood[node_index]][node_index][site] = 0.0;
                }
                else if ( use_missing_data == true && missing_data[left_index][site] == true && missing_data[right_index][site] == false )
                {
                    missing_data[node_index][site] = false;
                    
                    p_node[site]  = p_left[site] + p_right[site];
                    mu_node[site] = mu_right[site];
                    
                    this->contrast_uncertainty_per_site[this->active_likelihood[node_index]][node_index][site] = t_right;

                }
                else if ( use_missing_data == true && missing_data[left_index][site] == false && missing_data[right_index][site] == true )
                {
                    missing_data[node_index][site] = false;
                    
                    p_node[site]  = p_left[site] + p_right[site];
                    mu_node[site] = mu_left[site];
                    
                    this->contrast_uncertainty_per_site[this->active_likelihood[node_index]][node_index][site] = t_left;
                }
                else
                {

                    // get the site specific rate of evolution
                    double standDev = this->computeSiteRate(site) * stdev;

                    // compute the contrasts for this site and node
                    double contrast = mu_left[site] - mu_right[site];

                    // compute the probability for the contrasts at this node
                    double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);

                    // sum up the probabilities of the contrasts
                    p_node[site] = lnl_node + p_left[site] + p_right[site];
                    
                    mu_node[site] = (mu_left[site]*t_right + mu_right[site]*t_left) / (t_left+t_right);
                    
                    if ( use_missing_data == true )
                    {
                        missing_data[node_index][site] = false;
                        this->contrast_uncertainty_per_site[this->active_likelihood[node_index]][node_index][site] = (t_left*t_right) / (t_left+t_right);
                    }
                    
                }
                
            } // end for-loop over all sites
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}



void PhyloBrownianProcessStateDependent::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirty_nodes[index] )
    {
        // the root doesn't have an ancestor
        if ( !n.isRoot() )
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


void PhyloBrownianProcessStateDependent::redrawValue(void)
{
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new ContinuousCharacterData();
    
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxa = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    const Tree& tau = character_histories->getValue().getTree();
    
    // simulate the root sequence
    size_t root_index = tau.getRoot().getIndex();
    ContinuousTaxonData &root = taxa[ root_index ];
    
    std::vector<double> root_states = simulateRootCharacters(num_sites);
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // create the character
        double c = root_states[i];
        
        // add the character to the sequence
        root.addCharacter( c );
    }
    
    // recursively simulate the sequences
    simulateRecursively( tau.getRoot(), taxa );
    
    // we call now our method to resample the tips
    // this is important if we have multiple samples (e.g. individuals) per species
    simulateTipSamples( taxa );
    
    // tell the derived classes
    this->resetValue();

}



void PhyloBrownianProcessStateDependent::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods     = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrasts               = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    missing_data            = std::vector<std::vector<bool> >(this->num_nodes, std::vector<bool>(this->num_sites, false) );

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
    std::vector<TopologyNode*> nodes = this->character_histories->getValue().getTree().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                
                double c = taxon.getCharacter(site_indices[site]);
                
                if ( RbMath::isFinite(c) == false )
                {
                    missing_data[(*it)->getIndex()][site] = true;
                    use_missing_data = true;
                }
                
            }
        }
    }
    
    if ( use_missing_data == true )
    {
        contrast_uncertainty.clear();
        contrast_uncertainty_per_site   = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    }
    else
    {
        contrast_uncertainty_per_site.clear();
        contrast_uncertainty            = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );
    }
                
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                
                double c = taxon.getCharacter(site_indices[site]);
                
                if ( use_missing_data == false )
                {
                    contrast_uncertainty[0][(*it)->getIndex()] = 0;
                    contrast_uncertainty[1][(*it)->getIndex()] = 0;
                }
                else
                {
                    contrast_uncertainty_per_site[0][(*it)->getIndex()][site] = 0;
                    contrast_uncertainty_per_site[1][(*it)->getIndex()][site] = 0;
                }
                contrasts[0][(*it)->getIndex()][site] = c;
                contrasts[1][(*it)->getIndex()][site] = c;
            }
        }
    }
    
    
    // finally we set all the flags for recomputation
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = true;
    }
    
    // set the active likelihood pointers
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        active_likelihood[index] = 0;
        changed_nodes[index] = false;
    }
    
}



void PhyloBrownianProcessStateDependent::restoreSpecialization( const DagNode* affecter )
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


void PhyloBrownianProcessStateDependent::setSigma(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( state_dependent_sigma );
    homogeneous_sigma       = NULL;
    state_dependent_sigma   = NULL;
    
    
    // set the value
    homogeneous_sigma = s;
    
    // add the new parameter
    this->addParameter( homogeneous_sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloBrownianProcessStateDependent::setSigma(const TypedDagNode<RbVector<double> > *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( state_dependent_sigma );
    homogeneous_sigma       = NULL;
    state_dependent_sigma   = NULL;
    
    
    // set the value
    state_dependent_sigma = s;
    
    // add the new parameter
    this->addParameter( state_dependent_sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}



void PhyloBrownianProcessStateDependent::setValue(ContinuousCharacterData *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< ContinuousCharacterData >::setValue(v, force);
    
    // reset the number of sites
    num_sites = v->getNumberOfIncludedCharacters();
    
    // tell the derived classes
    this->resetValue();
    
}


void PhyloBrownianProcessStateDependent::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    
    // get the sequence of this node
    size_t node_index = node.getIndex();
    const ContinuousTaxonData &parent = taxa[ node_index ];
    
    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        
        // get the branch length for this child
        double branch_length = child.getBranchLength();
        
        // get the branch specific rate
        double branch_time = computeBranchTime( child.getIndex(), branch_length );
        
        //@TODO: Need to change the simulation to actually use the states from the character history
        
        // get the branch specific rate
        double branch_sigma = computeStateDependentSigma( 0 );
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double branch_rate = branch_time;
            double m = parent_state;
            
            double stand_dev = branch_sigma * sqrt(branch_rate);
            
            // create the character
            double c = RbStatistics::Normal::rv( m, stand_dev, *rng);
            
            // add the character to the sequence
            taxon.addCharacter( c );
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


std::vector<double> PhyloBrownianProcessStateDependent::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


void PhyloBrownianProcessStateDependent::simulateTipSamples( const std::vector< ContinuousTaxonData > &taxon_data )
{
    
    const Tree& tau = character_histories->getValue().getTree();
    // add the taxon data to the character data
    for (size_t i = 0; i < tau.getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxon_data[i] );
    }
    
}


double PhyloBrownianProcessStateDependent::sumRootLikelihood( void )
{
    // get the root node
    const Tree& tau = character_histories->getValue().getTree();
    const TopologyNode &root = tau.getRoot();
    
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


void PhyloBrownianProcessStateDependent::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    const Tree& tau = character_histories->getValue().getTree();

    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->state_dependent_sigma && false )
    {
        
        const std::set<size_t> &indices = this->state_dependent_sigma->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everyting for recomputation
            touchAll = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = tau.getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            }
        }
    }
    else // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
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
void PhyloBrownianProcessStateDependent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == homogeneous_sigma)
    {
        homogeneous_sigma = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == state_dependent_sigma)
    {
        state_dependent_sigma = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == character_histories)
    {
        character_histories = static_cast<const TypedDagNode< CharacterHistoryDiscrete >* >( newP );
    }
        
}


