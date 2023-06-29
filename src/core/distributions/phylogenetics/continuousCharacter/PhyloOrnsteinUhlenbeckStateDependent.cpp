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
#include "PhyloOrnsteinUhlenbeckStateDependent.h"
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

PhyloOrnsteinUhlenbeckStateDependent::PhyloOrnsteinUhlenbeckStateDependent(const TypedDagNode<Tree> *t, const StochasticNode<AbstractHomologousDiscreteCharacterData>* states, size_t ns) : AbstractPhyloContinuousCharacterProcess( t, ns ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrasts( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    contrast_uncertainty( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    normalizing_constants( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 1.0) ) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
    character_states( states )
{
    // initialize default parameters
    root_state                  = new ConstantNode<double>("", new double(0.0) );
    homogeneous_alpha           = new ConstantNode<double>("", new double(0.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_theta           = new ConstantNode<double>("", new double(0.0) );
    state_dependent_alpha       = NULL;
    state_dependent_sigma       = NULL;
    state_dependent_theta       = NULL;
    
    
    // add parameters
    addParameter( homogeneous_alpha );
    addParameter( homogeneous_sigma );
    addParameter( homogeneous_theta );
    addParameter( character_states );
    
    
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
PhyloOrnsteinUhlenbeckStateDependent::~PhyloOrnsteinUhlenbeckStateDependent( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloOrnsteinUhlenbeckStateDependent* PhyloOrnsteinUhlenbeckStateDependent::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckStateDependent( *this );
}


double PhyloOrnsteinUhlenbeckStateDependent::computeStateDependentAlpha(size_t state_idx) const
{
    
    // get the selection rate for the branch
    double a;
    if ( this->state_dependent_alpha != NULL )
    {
        a = this->state_dependent_alpha->getValue()[state_idx];
    }
    else
    {
        a = this->homogeneous_alpha->getValue();
    }
    
    return a;
}


double PhyloOrnsteinUhlenbeckStateDependent::computeStateDependentSigma(size_t state_idx) const
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


double PhyloOrnsteinUhlenbeckStateDependent::computeStateDependentTheta(size_t state_idx) const
{
    
    // get the optimum (theta) for the branch
    double t;
    if ( this->state_dependent_theta != NULL )
    {
        t = this->state_dependent_theta->getValue()[state_idx];
    }
    else
    {
        t = this->homogeneous_theta->getValue();
    }
    
    return t;
}


double PhyloOrnsteinUhlenbeckStateDependent::computeRootState( size_t state_index ) const
{
    
    // get the root-state parameter
    double root_state = this->root_state->getValue();
    
    return root_state;
}


double PhyloOrnsteinUhlenbeckStateDependent::computeLnProbability( void )
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
    {
                
        recursiveComputeLnProbability( root, root_index );
        
        // sum the partials up
        this->ln_prob = sumRootLikelihood();
        
    }
    return this->ln_prob;
}



void PhyloOrnsteinUhlenbeckStateDependent::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


void PhyloOrnsteinUhlenbeckStateDependent::keepSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckStateDependent::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    // check for recomputation
    if ( node.isTip() == false && dirty_nodes[node_index] == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        std::vector<double> &p_node             = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node            = this->contrasts[this->active_likelihood[node_index]][node_index];
        std::vector<double> &norm_const_node    = this->normalizing_constants[this->active_likelihood[node_index]][node_index];
        
        
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
            
            // get the per node and site normalizing constants
            const std::vector<double> &norm_const_left  = this->normalizing_constants[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &norm_const_right = this->normalizing_constants[this->active_likelihood[right_index]][right_index];
            
            // get the propagated uncertainties
            double delta_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
            double delta_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];
            
            
            // @TODO: maybe change later to allow for more states.
            double current_mu_left  = mu_left[0];
            double current_mu_right = mu_right[0];

            
            // get character history for tree
            const TreeHistoryCtmc<StandardState>* p = static_cast< const TreeHistoryCtmc<StandardState>* >(&character_states->getDistribution());
            
            // get branch history for node
            const BranchHistory& bh_left = p->getHistory(*left);
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_left = bh_left.getHistory();
            
//            const std::multiset<CharacterEvent*, CharacterEventCompare> &h_left = bh_left.getHistory();
//            CharacterEvent *event_left = *(h_left.begin());
//            size_t state_begin_left  = static_cast<CharacterEventDiscrete*>(bh_left.getChildCharacters()[0])->getState();
            
//            size_t previous_state_left = state_begin_left;
            double begin_time_left = left->getAge();
            
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_left.begin(); it!=hist_left.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                //            event_time = end_time;
                //            double time_interval = event_time - begin_time;
                
                // we need to set the current rate category
                size_t current_state = event->getState();
                
                
                double v_left  = 0;
                double time_left = event_time - begin_time_left;
                double sigma_left = computeStateDependentSigma( current_state );
                double alpha_left = computeStateDependentAlpha( current_state );
                double theta_left = computeStateDependentTheta( current_state );
                
                if ( alpha_left > 1E-20 )
                {
                    v_left = (sigma_left*sigma_left) / (2.0*alpha_left) * (exp(2.0*alpha_left*time_left) - 1.0 );
                }
                else
                {
                    v_left = (sigma_left*sigma_left) * time_left;
                }
                
                double var_left  = (v_left)  + delta_left  * exp(2.0*alpha_left *time_left);
                
                // now store the constrast variance for this part of the branch
                delta_left = var_left;
                
                double m_left   = exp(1.0 * time_left  * alpha_left ) * (current_mu_left  - theta_left)  + theta_left;
                current_mu_left = m_left;
                
                // update our variables
                begin_time_left     = event_time;
            }
            size_t last_state_left  = static_cast<CharacterEventDiscrete*>(bh_left.getParentCharacters()[0])->getState();
            
            
            
            
            // get branch history for node
            const BranchHistory& bh_right = p->getHistory(right);
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_right = bh_right.getHistory();
            
//            size_t state_begin_right  = static_cast<CharacterEventDiscrete*>(bh_right.getChildCharacters()[0])->getState();

//            size_t previous_state_right = state_begin_right;
            double begin_time_right = right.getAge();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_right.begin(); it!=hist_right.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                
                // we need to set the current rate category
                size_t current_state = event->getState();
                
                
                double v_right  = 0;
                double time_right = event_time - begin_time_right;
                double sigma_right = computeStateDependentSigma( current_state );
                double alpha_right = computeStateDependentAlpha( current_state );
                double theta_right = computeStateDependentTheta( current_state );
                
                if ( alpha_right > 1E-20 )
                {
                    v_right = (sigma_right*sigma_right) / (2.0*alpha_right) * (exp(2.0*alpha_right*time_right) - 1.0 );
                }
                else
                {
                    v_right = (sigma_right*sigma_right) * time_right;
                }
                
                double var_right  = (v_right)  + delta_right  * exp(2.0*alpha_right *time_right);
                
                // now store the constrast variance for this part of the branch
                delta_right = var_right;
                
                double m_right   = exp(1.0 * time_right  * alpha_right ) * (current_mu_right  - theta_right)  + theta_right;
                current_mu_right = m_right;
                
                // update our variables
                begin_time_right     = event_time;
//                previous_state_right = current_state;
            }
            size_t last_state_right  = static_cast<CharacterEventDiscrete*>(bh_right.getParentCharacters()[0])->getState();

            
            
            
            
            
            // get the scaled branch lengths
            double v_left  = 0;
//            if ( j == 1 )
//            {
                double bl_left = node.getAge() - begin_time_left;
                double sigma_left = computeStateDependentSigma(last_state_left);
                double alpha_left = computeStateDependentAlpha(last_state_left);
                if ( alpha_left > 1E-20 )
                {
                    v_left = (sigma_left*sigma_left) / (2.0*alpha_left) * (exp(2.0*alpha_left*bl_left) - 1.0 );
                }
                else
                {
                    v_left = (sigma_left*sigma_left) * bl_left;
                }
//            }
            
            double bl_right = node.getAge() - begin_time_right;
            double sigma_right = computeStateDependentSigma(last_state_right);
            double alpha_right = computeStateDependentAlpha(last_state_right);
            double v_right =0.0;
            if ( alpha_right > 1E-20 )
            {
                v_right = (sigma_right*sigma_right) / (2.0*alpha_right) * (exp(2.0*alpha_right*bl_right) - 1.0 );
            }
            else
            {
                v_right = (sigma_right*sigma_right) * bl_right;
            }
            
            // add the propagated uncertainty to the branch lengths
            double var_left  = (v_left)  + delta_left  * exp(2.0*alpha_left *bl_left);
            double var_right = (v_right) + delta_right * exp(2.0*alpha_right*bl_right);
            
            // set delta_node = (t_l*t_r)/(t_l+t_r);
            double var_node = (var_left*var_right) / (var_left+var_right);
            this->contrast_uncertainty[this->active_likelihood[node_index]][node_index] = var_node;
            
            double theta_left   = computeStateDependentTheta( last_state_left  );
            double theta_right  = computeStateDependentTheta( last_state_right );
            
            double stdev = sqrt(var_left+var_right);
            for (int i=0; i<this->num_sites; i++)
            {
                
                double m_left   = exp(1.0 * bl_left  * alpha_left ) * (current_mu_left  - theta_left)  + theta_left;
                double m_right  = exp(1.0 * bl_right * alpha_right) * (current_mu_right - theta_right) + theta_right;
                mu_node[i] = (m_left*var_right + m_right*var_left) / (var_left+var_right);
                
                // compute the contrasts for this site and node
                double contrast = m_left - m_right;
                
                // compute the probability for the contrasts at this node
                double z_left  = norm_const_left[i];
                double z_right = norm_const_right[i];
                double z_node  = exp(alpha_left*bl_left+alpha_right*bl_right)/(z_left*z_right) * exp( -1.0 * contrast * contrast / ( 2.0 *(var_left+var_right) ) ) / RbConstants::SQRT_2PI / stdev;
                norm_const_node[i] = 1.0;
                
                double lnl_node = log( z_node );
//                lnl_node -= RbConstants::LN_SQRT_2PI - log( sqrt( var_node ) );
//                lnl_node -= ( (contrast-mu_node[i])*(contrast-mu_node[i]) ) / (2.0*var_node);
//                lnl_node -= ( mu_node[i]*mu_node[i] ) / (2.0*var_node);
                
                if ( node.isRoot() == true )
                {
                    double root_state = computeRootState( last_state_left );
                    // dnorm(root.x, vals[1], sqrt(vals[2]), TRUE)
                    lnl_node += RbStatistics::Normal::lnPdf( root_state, sqrt((var_left*var_right) / (var_left+var_right)), mu_node[i]);
                }
                
                // sum up the probabilities of the contrasts
                p_node[i] = lnl_node + p_left[i] + p_right[i];
                
            } // end for-loop over all sites
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}



void PhyloOrnsteinUhlenbeckStateDependent::recursivelyFlagNodeDirty( const TopologyNode &n )
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


void PhyloOrnsteinUhlenbeckStateDependent::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods     = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrasts               = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrast_uncertainty    = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );
    normalizing_constants   = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 1.0) ) );

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
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(site_indices[site]);
                contrasts[0][(*it)->getIndex()][site] = c;
                contrasts[1][(*it)->getIndex()][site] = c;
                contrast_uncertainty[0][(*it)->getIndex()] = 0;
                contrast_uncertainty[1][(*it)->getIndex()] = 0;
                normalizing_constants[0][(*it)->getIndex()][site] = 1.0;
                normalizing_constants[1][(*it)->getIndex()][site] = 1.0;
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


void PhyloOrnsteinUhlenbeckStateDependent::restoreSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckStateDependent::setAlpha(const TypedDagNode<double> *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( state_dependent_alpha );
    homogeneous_alpha       = NULL;
    state_dependent_alpha   = NULL;
    
    
    // set the value
    homogeneous_alpha = a;
    
    // add the new parameter
    this->addParameter( homogeneous_alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckStateDependent::setAlpha(const TypedDagNode<RbVector<double> > *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( state_dependent_alpha );
    homogeneous_alpha       = NULL;
    state_dependent_alpha   = NULL;
    
    
    // set the value
    state_dependent_alpha   = a;
    
    // add the new parameter
    this->addParameter( state_dependent_alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckStateDependent::setRootState(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( root_state );
    root_state = s;
    
    // add the new parameter
    this->addParameter( root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckStateDependent::setSigma(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckStateDependent::setSigma(const TypedDagNode<RbVector<double> > *s)
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


void PhyloOrnsteinUhlenbeckStateDependent::setTheta(const TypedDagNode<double> *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( state_dependent_theta );
    homogeneous_theta        = NULL;
    state_dependent_theta    = NULL;
    
    
    // set the value
    homogeneous_theta = t;
    
    // add the new parameter
    this->addParameter( homogeneous_theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckStateDependent::setTheta(const TypedDagNode<RbVector<double> > *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( state_dependent_theta );
    homogeneous_theta        = NULL;
    state_dependent_theta    = NULL;
    
    
    // set the value
    state_dependent_theta = t;
    
    // add the new parameter
    this->addParameter( state_dependent_theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckStateDependent::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        
        // get the branch specific optimum (theta)
        double branch_theta = computeStateDependentTheta( 0 );
        
        // get the branch specific optimum (theta)
        double branch_alpha = computeStateDependentAlpha( 0 );
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double branch_rate = branch_time;
            
            double e = exp(-branch_alpha * branch_rate);
            double e2 = exp(-2.0 * branch_alpha * branch_rate);
            double m = e * parent_state + (1 - e) * branch_theta;
            
            double stand_dev = 0.0;
            if ( branch_alpha > 1E-10 )
            {
                double sigma_square = branch_sigma * branch_sigma;
                stand_dev = sqrt( (sigma_square / (2.0*branch_alpha)*(1.0 - e2)) );
            }
            else
            {
                // compute the standard deviation for this site
                stand_dev = branch_sigma * sqrt(branch_rate);
            }
            
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


std::vector<double> PhyloOrnsteinUhlenbeckStateDependent::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        // @TODO: need to take the root states from the character history
        chars[i] = computeRootState( 0 );
    }
    
    return chars;
}


double PhyloOrnsteinUhlenbeckStateDependent::sumRootLikelihood( void )
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


void PhyloOrnsteinUhlenbeckStateDependent::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
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
void PhyloOrnsteinUhlenbeckStateDependent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    
    if (oldP == root_state)
    {
        root_state = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == homogeneous_alpha)
    {
        homogeneous_alpha = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == state_dependent_alpha)
    {
        state_dependent_alpha = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_sigma)
    {
        homogeneous_sigma = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == state_dependent_sigma)
    {
        state_dependent_sigma = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_theta)
    {
        homogeneous_theta = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == state_dependent_theta)
    {
        state_dependent_theta = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == character_states)
    {
        character_states = static_cast<const StochasticNode< AbstractHomologousDiscreteCharacterData >* >( newP );
    }
    
    this->AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    
}


