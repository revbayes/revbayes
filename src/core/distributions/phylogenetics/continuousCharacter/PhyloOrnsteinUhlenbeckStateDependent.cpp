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

PhyloOrnsteinUhlenbeckStateDependent::PhyloOrnsteinUhlenbeckStateDependent(const TypedDagNode<CharacterHistoryDiscrete> *ch, size_t ns) : TypedDistribution< ContinuousCharacterData > ( new ContinuousCharacterData() ),
    num_nodes( ch->getValue().getNumberBranches()+1 ),
    num_sites( ns ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    means( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    variances( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
    character_histories( ch )
{
    // initialize default parameters
    root_state                  = new ConstantNode<double>("", new double(0.0) );
    homogeneous_alpha           = new ConstantNode<double>("", new double(0.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_theta           = new ConstantNode<double>("", new double(0.0) );
    state_dependent_alpha       = NULL;
    state_dependent_sigma       = NULL;
    state_dependent_theta       = NULL;
    observation_variance        = new ConstantNode<double>("", new double(0.0) );
    
    
    // add parameters
    addParameter( homogeneous_alpha );
    addParameter( homogeneous_sigma );
    addParameter( homogeneous_theta );
    addParameter( observation_variance );
    addParameter( character_histories );
    addParameter( root_state );
    
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
    
//    // remove myself from the tree listeners
//    if ( tau != NULL )
//    {
//        tau->getValue().getTreeChangeEventHandler().removeListener( this );
//    }
    
}



PhyloOrnsteinUhlenbeckStateDependent* PhyloOrnsteinUhlenbeckStateDependent::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckStateDependent( *this );
}



double PhyloOrnsteinUhlenbeckStateDependent::computeBranchTime( size_t nide_idx, double brlen )
{
    
    // get the clock rate for the branch
    double branch_time = 1.0;
//    if ( this->heterogeneous_clock_rates != NULL )
//    {
//        double sigma = this->heterogeneous_clock_rates->getValue()[nide_idx];
//        branch_time = sigma * sigma * brlen;
//    }
//    else
//    {
//        double sigma = this->homogeneous_clock_rate->getValue();
//        branch_time = sigma * sigma * brlen;
//    }
    
    return branch_time;
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

double PhyloOrnsteinUhlenbeckStateDependent::computeObservationVariance() const
{
    double ov = this->observation_variance->getValue();
    return ov;
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
    double root_state_value = this->root_state->getValue();
    
    return root_state_value;
}


double PhyloOrnsteinUhlenbeckStateDependent::computeLnProbability( void )
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

    std::vector<TopologyNode*> nodes = tau.getNodes();

    double ov;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        ov = computeObservationVariance();

        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                //variances[0][(*it)->getIndex()] = ov;
                //variances[1][(*it)->getIndex()] = ov;
                size_t node_index = (*it)->getIndex();
                this->variances[this->active_likelihood[node_index]][node_index] = ov;
            }
        }
    }       

    // only necessary if the root is actually dirty
    //if ( this->dirty_nodes[root_index] )
    //{
                
        recursiveComputeLnProbability( root, root_index );
        
        // sum the partials up
        this->ln_prob = sumRootLikelihood();
        
    //    }
   
    //std::cout << this->ln_prob << std::endl;
    //std::cout << "observation variance: \t" << ov << std::endl;
    //std::cout << "log likelihood: \t" << this->ln_prob << std::endl;
    //std::cout.flush();
    return this->ln_prob;
}



//void PhyloOrnsteinUhlenbeckStateDependent::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
//{
//
//    // call a recursive flagging of all node above (closer to the root) and including this node
//    recursivelyFlagNodeDirty( n );
//
//}


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
        
        std::vector<double> &p_node       = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node                   = this->means[this->active_likelihood[node_index]][node_index];

        
        
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
            
            // get the means for the left and right subtrees
            const std::vector<double> &mu_left  = this->means[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->means[this->active_likelihood[right_index]][right_index];
            
            // get the variances of the left and right child nodes
            double delta_left  = this->variances[this->active_likelihood[left_index]][left_index];
            double delta_right = this->variances[this->active_likelihood[right_index]][right_index];
            
            
            // @TODO: maybe change later to allow for more characters/sites.
            double mean_left  = mu_left[0];
            double mean_right = mu_right[0];
            
            // get character history for tree
            const CharacterHistory& current_history = character_histories->getValue();
            
            // get branch history for the left branch
            const BranchHistory& bh_left = current_history.getHistory(left_index);
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_left = bh_left.getHistory();
            
            double v_left;
            double log_nf_left = 0;
            double begin_time_left = left->getAge();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_left.begin(); it!=hist_left.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                // we need to find the current state
                size_t current_state = event->getState();
                
                double time_left = event_time - begin_time_left;
                double sigma_left = computeStateDependentSigma( current_state );
                double alpha_left = computeStateDependentAlpha( current_state );
                double theta_left = computeStateDependentTheta( current_state );
               
                if ( alpha_left > 1E-20 )
                {
                    v_left  = (sigma_left*sigma_left) / (2.0*alpha_left) * (exp(2.0*alpha_left*time_left) - 1.0 );
                    mean_left  = exp(1.0 * time_left  * alpha_left ) * (mean_left  - theta_left)  + theta_left;
                }
                else
                {
                    v_left  = (sigma_left*sigma_left) * time_left;
                    //mean_left  = current_mu_left; // no change
                }
                delta_left = v_left + delta_left * exp(2.0*alpha_left *time_left);
                begin_time_left = event_time;
                
                // update the log normalizing factor
                log_nf_left += time_left * alpha_left;
            }
            size_t last_state_left  = static_cast<CharacterEventDiscrete*>(bh_left.getParentCharacters()[0])->getState();
            
            
            // get branch history for the right branch
            const BranchHistory& bh_right = current_history.getHistory(right_index);
            const std::multiset<CharacterEvent*,CharacterEventCompare>& hist_right = bh_right.getHistory();

            double v_right;
            double log_nf_right = 0;
            double begin_time_right = right.getAge();
            for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_iterator it=hist_right.begin(); it!=hist_right.end(); ++it)
            {
                CharacterEventDiscrete* event = static_cast<CharacterEventDiscrete*>(*it);
                double event_time = event->getAge();
                
                // we need to set the current rate category
                size_t current_state = event->getState();
                
                double time_right = event_time - begin_time_right;
                double sigma_right = computeStateDependentSigma( current_state );
                double alpha_right = computeStateDependentAlpha( current_state );
                double theta_right = computeStateDependentTheta( current_state );
               
                if ( alpha_right > 1E-20 )
                {
                    v_right = (sigma_right*sigma_right) / (2.0*alpha_right) * (exp(2.0*alpha_right*time_right) - 1.0 );
                    mean_right   = exp(1.0 * time_right  * alpha_right ) * (mean_right  - theta_right)  + theta_right;
                }
                else
                {
                    v_right = (sigma_right*sigma_right) * time_right;
                    //mean_right = current_mu_right; no change
                }
                delta_right = v_right + delta_right * exp(2.0*alpha_right *time_right);
               
                // update the time
                begin_time_right     = event_time;
                // update the log normalizing factor
                log_nf_right += time_right * alpha_right;
            }
            size_t last_state_right  = static_cast<CharacterEventDiscrete*>(bh_right.getParentCharacters()[0])->getState();


            // the above code does not run if there are no transitions. In either case, we need to finish the "final" branch segment
            double bl_left = node.getAge() - begin_time_left;
            double sigma_left = computeStateDependentSigma(last_state_left);
            double alpha_left = computeStateDependentAlpha(last_state_left);
            double theta_left = computeStateDependentTheta(last_state_left);
            if ( alpha_left > 1E-20 )
            {
                v_left = (sigma_left*sigma_left) / (2.0*alpha_left) * (exp(2.0*alpha_left*bl_left) - 1.0 );
                mean_left  = exp(1.0 * bl_left  * alpha_left ) * (mean_left  - theta_left)  + theta_left;
            }
            else
            {
                v_left = (sigma_left*sigma_left) * bl_left;
                //mean_left = nothing to change
            }
            delta_left = v_left + delta_left * exp(2.0*alpha_left *bl_left);
            log_nf_left += bl_left * alpha_left;
            
            // same but for right branch
            double bl_right = node.getAge() - begin_time_right;
            double sigma_right = computeStateDependentSigma(last_state_right);
            double alpha_right = computeStateDependentAlpha(last_state_right);
            double theta_right = computeStateDependentTheta(last_state_right);
            if ( alpha_right > 1E-20 )
            {
                v_right = (sigma_right*sigma_right) / (2.0*alpha_right) * (exp(2.0*alpha_right*bl_right) - 1.0 );
                mean_right  = exp(1.0 * bl_right  * alpha_right ) * (mean_right  - theta_right)  + theta_right;
            }
            else
            {
                v_right = (sigma_right*sigma_right) * bl_right;
                //mean_right = nothing to change
            }
            delta_right = v_right + delta_right * exp(2.0*alpha_right *bl_right);
            log_nf_right += bl_right * alpha_right;

            double var_left = delta_left;
            double var_right = delta_right;
            
            // calculate and store the node variance
            double var_node = (var_left*var_right) / (var_left+var_right);
            this->variances[this->active_likelihood[node_index]][node_index] = var_node;


            // this loop is broken currently, does not work for multiple characters
            for (int i=0; i<this->num_sites; i++)
            {
                mu_node[i] = (var_left*mean_right + var_right*mean_left) / (var_left+var_right);
                
                // compute the contrasts for this site and node
                double contrast = mean_left - mean_right;
                
                double a = -(contrast*contrast / (2*(var_left + var_right)));               
                double b = log(2*RbConstants::PI*(var_left+var_right))/2.0;
                double lnl_node = log_nf_left + log_nf_right + a - b;
                
                if ( node.isRoot() == true )
                {
                    double root_state = computeRootState( last_state_left );
                    lnl_node += RbStatistics::Normal::lnPdf( root_state, sqrt(var_node), mu_node[i]);
                }
                // sum up the log normalizing factors of the subtrees
                p_node[i] = lnl_node + p_left[i] + p_right[i];
                
            } // end for-loop over all sites
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.

    //if ( node.isTip() == true && dirty_nodes[node_index] == true ){
     //   double ov = computeObservationVariance(); 
      //  this->variances[this->active_likelihood[node_index]][node_index] = ov;
       // dirty_nodes[node_index] = false; // why need to do this?
        
        //std::cout << "observation variance: \t" << ov << std::endl;
        //std::cout << "sigma state 1: \t" << computeStateDependentSigma(0) << std::endl;
        //std::cout.flush();
//    }
    
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


void PhyloOrnsteinUhlenbeckStateDependent::redrawValue(void)
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



void PhyloOrnsteinUhlenbeckStateDependent::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods  = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    means                = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    variances            = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );


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
    
    const Tree& tau = character_histories->getValue().getTree();
    std::vector<TopologyNode*> nodes = tau.getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        double ov = computeObservationVariance();

        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(site_indices[site]);
                means[0][(*it)->getIndex()][site] = c;
                means[1][(*it)->getIndex()][site] = c;
                variances[0][(*it)->getIndex()] = ov;
                variances[1][(*it)->getIndex()] = ov;
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

void PhyloOrnsteinUhlenbeckStateDependent::setObservationVariance(const TypedDagNode<double> *ov)
{
    // remove the old parameter first
    this->removeParameter( observation_variance );
    
    // set the value
    observation_variance = ov;
    
    // add the new parameter
    this->addParameter( observation_variance );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}

void PhyloOrnsteinUhlenbeckStateDependent::setValue(ContinuousCharacterData *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< ContinuousCharacterData >::setValue(v, force);
    
    // reset the number of sites
    num_sites = v->getNumberOfIncludedCharacters();
    
    // tell the derived classes
    this->resetValue();
    
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


void PhyloOrnsteinUhlenbeckStateDependent::simulateTipSamples( const std::vector< ContinuousTaxonData > &taxon_data )
{
    
    const Tree& tau = character_histories->getValue().getTree();
    // add the taxon data to the character data
    for (size_t i = 0; i < tau.getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxon_data[i] );
    }
    
}


double PhyloOrnsteinUhlenbeckStateDependent::sumRootLikelihood( void )
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


void PhyloOrnsteinUhlenbeckStateDependent::touchSpecialization( const DagNode* affecter, bool touchAll )
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
void PhyloOrnsteinUhlenbeckStateDependent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == root_state)
    {
        root_state = static_cast<const TypedDagNode< double >* >( newP );
    }

    if (oldP == observation_variance)
    {
        observation_variance = static_cast<const TypedDagNode< double >* >( newP );
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
    
    if (oldP == character_histories)
    {
        character_histories = static_cast<const TypedDagNode< CharacterHistoryDiscrete >* >( newP );
    }
        
}


