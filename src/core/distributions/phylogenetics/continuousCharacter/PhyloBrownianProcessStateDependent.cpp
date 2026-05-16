#include <cmath>
#include <cstddef>
#include <deque>
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
    root_value                  = new ConstantNode<double>("", new double(0.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    state_dependent_sigma       = NULL;


    // add parameters
    addParameter( homogeneous_sigma );
    addParameter( character_histories );
    addParameter( root_value );

    // now we need to reset the value
    this->redrawValue();

    // we need to reset the contrasts
    resetValue();
}


/**
 * Destructor.
 */
PhyloBrownianProcessStateDependent::~PhyloBrownianProcessStateDependent( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

}



PhyloBrownianProcessStateDependent* PhyloBrownianProcessStateDependent::clone( void ) const
{

    return new PhyloBrownianProcessStateDependent( *this );
}


double PhyloBrownianProcessStateDependent::computeRootValue( void ) const
{
    // get the root-value parameter
    double rvl = this->root_value->getValue();

    return rvl;
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


LogDensity PhyloBrownianProcessStateDependent::computeLnProbability( void )
{

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

                    // compute the contrasts for this site and node
                    double contrast = mu_left[site] - mu_right[site];

                    // compute the probability for the contrasts at this node
                    double lnl_node = RbStatistics::Normal::lnPdf(0, stdev, contrast);

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
                throw RbException( "The number of sites specified does not match the number of characters included." );
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


void PhyloBrownianProcessStateDependent::setRootValue(const TypedDagNode<double> *rvl)
{

    // remove the old parameter first
    this->removeParameter( root_value );
    root_value = NULL;

    // set the value
    root_value = rvl;

    // add the new parameter
    this->addParameter( root_value );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
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

    size_t number_states = character_histories->getValue().getNumberOfStates();
    // make sure that the state-space is correct
    if ( s->getValue().size() != number_states )
    {
        throw RbException() << "The number of states (" << number_states << ") in the character history doesn't match the number of sigma parameters (" << s->getValue().size() << ")";
    }


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


double PhyloBrownianProcessStateDependent::simulateEpisode(size_t state_index, double delta_t, double ancestral_value)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the parameter values
    double sigma = computeStateDependentSigma(state_index);

    // calculate the mean and the variance
    double mu = ancestral_value;
    double sd = sigma * sqrt(delta_t);

    // draw the new character state as a Gaussian random variable
    double y = RbStatistics::Normal::rv(mu, sd, *rng);

    return y;
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

        // get the branch history
        const CharacterHistory& current_history = character_histories->getValue();
        size_t child_index = child.getIndex();
        const BranchHistory& bh = current_history.getHistory(child_index);

        const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh.getHistory();

        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            double youngest_time = child.getAge();
            double begin_time = youngest_time;

            // the episode states and times (if there was at least one discrete character state change)
            // the loop is from young to old
            // since it's pushed to the front of the deque,
            // the array is in order of old to young

            std::deque<double> times;
            std::deque<size_t> states;

            for (std::multiset<CharacterEvent*, CharacterEventCompare>::const_iterator iter = history.begin(); iter != history.end(); ++iter)
            {
                // get the state change event
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*iter);

                // calculate the times
                double event_time = event->getAge();
                double delta_t = event_time - begin_time;
                begin_time = event_time;

                // get the state index
                size_t current_state = event->getState();

                // save the
                times.push_front(delta_t);
                states.push_front(current_state);
            }

            // do it again, since the iterator above only does n-1 of the episodes
            size_t first_state = static_cast<CharacterEventDiscrete*>(bh.getParentCharacters()[0])->getState();
            double first_delta_t = node.getAge() - begin_time;

            times.push_front(first_delta_t);
            states.push_front(first_state);

            // get the ancestral character for this site
            double y = parent.getCharacter( i );

            // simulate the episodes
            for (size_t j = 0; j < times.size(); ++j)
            {
                size_t state = states[j];
                double delta_t = times[j];
                y = simulateEpisode(state, delta_t, y);
            }

            taxon.addCharacter(y);
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
        chars[i] = computeRootValue();
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

    // currently we don't actually use any fetching of parts of the tree that don't need recomputation.
    // this code is merely as a placeholder when this implemented in the future here.
    bool touch_all = true;

    if ( affecter == this->dag_node )
    {
        resetValue();
    }


    if ( touch_all )
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
