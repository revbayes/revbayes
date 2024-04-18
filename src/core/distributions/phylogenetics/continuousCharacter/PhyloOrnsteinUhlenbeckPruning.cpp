#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "PhyloOrnsteinUhlenbeckPruning.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloContinuousCharacterProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloOrnsteinUhlenbeckPruning::PhyloOrnsteinUhlenbeckPruning(const TypedDagNode<Tree> *t, size_t ns) : AbstractPhyloContinuousCharacterProcess( t, ns ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    means( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    variances( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) )
{
    // initialize default parameters
    root_state                  = new ConstantNode<double>("", new double(0.0) );
    homogeneous_alpha           = new ConstantNode<double>("", new double(0.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_theta           = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_alpha         = NULL;
    heterogeneous_sigma         = NULL;
    heterogeneous_theta         = NULL;
    
    
    // add parameters
    addParameter( homogeneous_alpha );
    addParameter( homogeneous_sigma );
    addParameter( homogeneous_theta );
    
    
    // We don'e want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // now we need to reset the value
    this->redrawValue();
    
    // we need to reset the means and variances
    resetValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloOrnsteinUhlenbeckPruning::~PhyloOrnsteinUhlenbeckPruning( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloOrnsteinUhlenbeckPruning* PhyloOrnsteinUhlenbeckPruning::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckPruning( *this );
}


double PhyloOrnsteinUhlenbeckPruning::computeBranchAlpha(size_t branch_idx) const
{
    
    // get the selection rate for the branch
    double a;
    if ( this->heterogeneous_alpha != NULL )
    {
        a = this->heterogeneous_alpha->getValue()[branch_idx];
    }
    else
    {
        a = this->homogeneous_alpha->getValue();
    }
    
    return a;
}


double PhyloOrnsteinUhlenbeckPruning::computeBranchSigma(size_t branch_idx) const
{
    
    // get the drift rate for the branch
    double s;
    if ( this->heterogeneous_sigma != NULL )
    {
        s = this->heterogeneous_sigma->getValue()[branch_idx];
    }
    else
    {
        s = this->homogeneous_sigma->getValue();
    }
    
    return s;
}


double PhyloOrnsteinUhlenbeckPruning::computeBranchTheta(size_t branch_idx) const
{
    
    // get the optimum (theta) for the branch
    double t;
    if ( this->heterogeneous_theta != NULL )
    {
        t = this->heterogeneous_theta->getValue()[branch_idx];
    }
    else
    {
        t = this->homogeneous_theta->getValue();
    }
    
    return t;
}


double PhyloOrnsteinUhlenbeckPruning::computeRootState( void ) const
{
    
    // get the root-state parameter
    double root_state = this->root_state->getValue();
    
    return root_state;
}


double PhyloOrnsteinUhlenbeckPruning::computeLnProbability( void )
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
    size_t rootIndex = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( this->dirty_nodes[rootIndex] )
    {
        
        
        recursiveComputeLnProbability( root, rootIndex );
        
        // sum the partials up
        this->ln_prob = sumRootLikelihood();
        
    }
    return this->ln_prob;
}



void PhyloOrnsteinUhlenbeckPruning::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


void PhyloOrnsteinUhlenbeckPruning::keepSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckPruning::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    // check for recomputation
    if ( node.isTip() == false && dirty_nodes[node_index] == true )
    {
        // mark as computed
        dirty_nodes[node_index] = false;
        
        std::vector<double> &p_node             = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double> &mu_node            = this->means[this->active_likelihood[node_index]][node_index];
        
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
            
            // get the means for the left and right subtrees
            const std::vector<double> &mu_left  = this->means[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->means[this->active_likelihood[right_index]][right_index];
            
            // get the variances for the left and right subtrees
            // the name "delta" comes from Felsenstein (2004; Fig.23.3, page 407)
            double delta_left  = this->variances[this->active_likelihood[left_index]][left_index];
            double delta_right = this->variances[this->active_likelihood[right_index]][right_index];
            
            // calculate the variance accounting for the branch
            double v_left  = 0.0;
            double bl_left = left->getBranchLength();
            double sigma_left = computeBranchSigma(left_index);
            double alpha_left = computeBranchAlpha(left_index);
            if ( alpha_left > 1E-20 )
            {
                v_left = (sigma_left*sigma_left) / (2.0*alpha_left) * expm1(2.0*alpha_left*bl_left);
            }
            else
            {
                v_left = (sigma_left*sigma_left) * bl_left;
            }
        
            double v_right = 0.0;
            double bl_right = right.getBranchLength();
            double sigma_right = computeBranchSigma(right_index);
            double alpha_right = computeBranchAlpha(right_index);
            if ( alpha_right > 1E-20 )
            {
                v_right = (sigma_right*sigma_right) / (2.0*alpha_right) * expm1(2.0*alpha_right*bl_right);
            }
            else
            {
                v_right = (sigma_right*sigma_right) * bl_right;
            }
            
            // add the branch variance with the variance of the subtrees
            //                 (this branch)       (                 subtree                )
            double var_left  = v_left            + delta_left  * exp(2.0*alpha_left *bl_left);
            double var_right = v_right           + delta_right * exp(2.0*alpha_right*bl_right);
            
            double theta_left   = computeBranchTheta( left_index );
            double theta_right  = computeBranchTheta( right_index );
           
            // the next steps are setting up the Gaussian variable
            // according to the steps outlined in the supplement of 
            // FitzJohn (2012, Methods in Ecol Evol). The Gaussian variable
            // has three components (equation 6)
            //
            // i) a variance
            // ii) a mean
            // iii) a normalizing factor

            
            // i) computing the variance at the node
            double var_node = (var_left*var_right) / (var_left+var_right);
            this->variances[this->active_likelihood[node_index]][node_index] = var_node;
            
            for (int i=0; i<this->num_sites; i++)
            {
                // ii) computing the mean at the node
                double mean_left   = exp(1.0 * bl_left  * alpha_left ) * (mu_left[i]  - theta_left)  + theta_left;
                double mean_right  = exp(1.0 * bl_right * alpha_right) * (mu_right[i] - theta_right) + theta_right;
                double mean_node = (mean_left*var_right + mean_right*var_left) / (var_left+var_right);
                mu_node[i] = mean_node; 
                
                
                // iii) computing the (log) normalizing factor at the node
                // 
                // unlike in FitzJohn(2012) supplement eq. (6), we don't propagate the 
                // normalizing factors (1/z_left) and (1/z_right) of the left and
                // right subtrees, since this can lead to underflow. Instead we store the log
                // normalizing factors, and add them to the partial log likelihoods
                double log_nf_left  = alpha_left  * bl_left;
                double log_nf_right = alpha_right * bl_right;
                double contrast = mean_left - mean_right;
                double a = -1.0 * contrast * contrast / ( 2.0 *(var_left+var_right) );
                double b = 0.5 * log( 2*RbConstants::PI*(var_left+var_right) );

                double log_norm_factor = log_nf_left + log_nf_right + a - b;
                double lnl_node = log_norm_factor;
                
                if ( node.isRoot() == true )
                {
                    // if the assumed root state is equal to the optimum at the root, then
                    // this pruning algorithm is 100% equivalent to the likelihood obtained
                    // using generalized least squares with a variance-covariance matrix
                    // for the residuals (r = y - theta), also called the vcv-method (introduced in Hansen 1997)
                    double root_state = computeRootState();
                    lnl_node += RbStatistics::Normal::lnPdf( root_state, sqrt(var_node), mu_node[i]);
                }
                
                // sum up the probabilities of the subtrees
                p_node[i] = lnl_node + p_left[i] + p_right[i];
                
            } // end for-loop over all sites
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}



void PhyloOrnsteinUhlenbeckPruning::recursivelyFlagNodeDirty( const TopologyNode &n )
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


void PhyloOrnsteinUhlenbeckPruning::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    means = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    variances = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );

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
                means[0][(*it)->getIndex()][site] = c;
                means[1][(*it)->getIndex()][site] = c;
                variances[0][(*it)->getIndex()] = 0;
                variances[1][(*it)->getIndex()] = 0;
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


void PhyloOrnsteinUhlenbeckPruning::restoreSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckPruning::setAlpha(const TypedDagNode<double> *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( heterogeneous_alpha );
    homogeneous_alpha      = NULL;
    heterogeneous_alpha    = NULL;
    
    
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


void PhyloOrnsteinUhlenbeckPruning::setAlpha(const TypedDagNode<RbVector<double> > *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( heterogeneous_alpha );
    homogeneous_alpha      = NULL;
    heterogeneous_alpha    = NULL;
    
    
    // set the value
    heterogeneous_alpha = a;
    
    // add the new parameter
    this->addParameter( heterogeneous_alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckPruning::setRootState(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckPruning::setSigma(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( heterogeneous_sigma );
    homogeneous_sigma      = NULL;
    heterogeneous_sigma    = NULL;
    
    
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


void PhyloOrnsteinUhlenbeckPruning::setSigma(const TypedDagNode<RbVector<double> > *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( heterogeneous_sigma );
    homogeneous_sigma      = NULL;
    heterogeneous_sigma    = NULL;
    
    
    // set the value
    heterogeneous_sigma = s;
    
    // add the new parameter
    this->addParameter( heterogeneous_sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckPruning::setTheta(const TypedDagNode<double> *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( heterogeneous_theta );
    homogeneous_theta      = NULL;
    heterogeneous_theta    = NULL;
    
    
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


void PhyloOrnsteinUhlenbeckPruning::setTheta(const TypedDagNode<RbVector<double> > *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( heterogeneous_theta );
    homogeneous_theta      = NULL;
    heterogeneous_theta    = NULL;
    
    
    // set the value
    heterogeneous_theta = t;
    
    // add the new parameter
    this->addParameter( heterogeneous_theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckPruning::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        
        // get the branch specific rate
        double branch_sigma = computeBranchSigma( child.getIndex() );
        
        // get the branch specific optimum (theta)
        double branch_theta = computeBranchTheta( child.getIndex() );
        
        // get the branch specific optimum (theta)
        double branch_alpha = computeBranchAlpha( child.getIndex() );
        
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


std::vector<double> PhyloOrnsteinUhlenbeckPruning::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState();
    }
    
    return chars;
}


double PhyloOrnsteinUhlenbeckPruning::sumRootLikelihood( void )
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


void PhyloOrnsteinUhlenbeckPruning::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_sigma )
    {
        
        const std::set<size_t> &indices = this->heterogeneous_sigma->getTouchedElementIndices();
        
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
void PhyloOrnsteinUhlenbeckPruning::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
    else if (oldP == heterogeneous_alpha)
    {
        heterogeneous_alpha = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_sigma)
    {
        homogeneous_sigma = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_sigma)
    {
        heterogeneous_sigma = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_theta)
    {
        homogeneous_theta = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_theta)
    {
        heterogeneous_theta = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    this->AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    
}


