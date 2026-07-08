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
#include "TreeChangeEventMessage.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloOrnsteinUhlenbeckPruning::PhyloOrnsteinUhlenbeckPruning(const TypedDagNode<Tree> *t, size_t ns) : AbstractPhyloContinuousCharacterProcess( t, ns ),
    node_likelihoods(this->num_nodes)
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
        resetValue();
    }
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();
    
    // only necessary if the root cache is invalid
    if ( not this->node_likelihoods.is_valid(rootIndex) )
    {
        
        
        recursiveComputeLnProbability( root, rootIndex );
        
        // sum the partials up
        this->ln_prob = sumRootLikelihood();
        
    }
    return this->ln_prob;
}



void PhyloOrnsteinUhlenbeckPruning::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    if (m == TreeChangeEventMessage::BRANCH_LENGTH)
    {
        invalidateBranchAndAncestors( n );
    }
    else
    {
        invalidateInternalNodes();
    }
}

// this function changes mu, variance and log_nf in-place
//
// calculate the variance accounting for the branch
// the next steps are setting up the Gaussian variable
// according to the steps outlined in the supplement of
// FitzJohn (2012, Methods in Ecol Evol). The Gaussian variable
// has three components (equation 6)
//
// i) a mean
// ii) a variance
// iii) a normalizing factor
void PhyloOrnsteinUhlenbeckPruning::propagateAuxiliaryVariables(double &mu, double &variance, double &log_nf, const TopologyNode& node )
{
    size_t node_index = node.getIndex();
    double time = node.getBranchLength();

    double theta = computeBranchTheta(node_index);
    double sigma = computeBranchSigma(node_index);
    double alpha = computeBranchAlpha(node_index);
               
    double v;
    if ( alpha > 1E-20 )
    {
        v = (sigma*sigma) / (2.0*alpha) * (exp(2.0*alpha*time) - 1.0 );
        mu  = exp(1.0 * time  * alpha ) * (mu  - theta)  + theta;
    }
    else
    {
        v  = (sigma*sigma) * time;
    }
    variance = v + variance * exp(2.0*alpha *time);
                
    // update the log normalizing factor
    log_nf += time * alpha;
}


void PhyloOrnsteinUhlenbeckPruning::keepSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.keep();
    }
}


void PhyloOrnsteinUhlenbeckPruning::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{

    // check for recomputation
    if ( node.isTip() == false && not node_likelihoods.is_valid(node_index) )
    {
        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        if (num_children != 2 )
        {
            throw RbException("internal node in the phylogeny does not have two descendants (in PhyloOrnsteinUhlenbeckPruning), not supported");
        }
            
        const TopologyNode &left = node.getChild(0);
        size_t left_index = left.getIndex();
        recursiveComputeLnProbability( left, left_index );

        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        recursiveComputeLnProbability( right, right_index );

        NodeLikelihood &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.means.resize(this->num_sites);
        node_cache.variances.resize(this->num_sites);
        node_cache.partial_likelihoods.resize(this->num_sites);
        node_cache.missing_data.resize(this->num_sites);

        std::vector<double> &mu_node            = node_cache.means;
        std::vector<double> &v_node             = node_cache.variances;
        std::vector<double> &p_node             = node_cache.partial_likelihoods;
        std::vector<bool> &missing_node         = node_cache.missing_data;

        // get the means for the left and right subtrees
        const NodeLikelihood &left_cache  = this->node_likelihoods[left_index];
        const NodeLikelihood &right_cache = this->node_likelihoods[right_index];
        const std::vector<double> &mu_left  = left_cache.means;
        const std::vector<double> &mu_right = right_cache.means;

        const std::vector<double> &v_left   = left_cache.variances;
        const std::vector<double> &v_right  = right_cache.variances;

        const std::vector<double> &p_left   = left_cache.partial_likelihoods;
        const std::vector<double> &p_right  = right_cache.partial_likelihoods;
        const std::vector<bool> &missing_left  = left_cache.missing_data;
        const std::vector<bool> &missing_right = right_cache.missing_data;
        
        
        size_t num_sites = this->num_sites;

        for (size_t i=0; i < num_sites; i++)
        {
            bool left_missing = missing_left[i];
            bool right_missing = missing_right[i];

            if ( use_missing_data == true && left_missing && right_missing )
            {
                missing_node[i] = true;

                mu_node[i] = RbConstants::Double::nan;
                v_node[i]  = 0.0;
                p_node[i]  = p_left[i] + p_right[i];
            }
            else if ( use_missing_data == true && left_missing && !right_missing )
            {
                missing_node[i] = false;
                
                double mean_right = mu_right[i];
                double var_right  = v_right[i];
                double log_nf_right = p_right[i];

                propagateAuxiliaryVariables(mean_right, var_right, log_nf_right, right);

                mu_node[i] = mean_right;
                v_node[i] = var_right;
                p_node[i] = log_nf_right;

            }
            else if ( use_missing_data == true && !left_missing && right_missing )
            {
                missing_node[i] = false;
                
                double mean_left = mu_left[i];
                double var_left = v_left[i];
                double log_nf_left = p_left[i];

                propagateAuxiliaryVariables(mean_left, var_left, log_nf_left, left);

                mu_node[i] = mean_left;
                v_node[i] = var_left;
                p_node[i] = log_nf_left;

            }
            else
            {
                // update mu, v and z for the left branch
                double mean_left = mu_left[i];
                double log_nf_left = p_left[i];
                double var_left = v_left[i];

                propagateAuxiliaryVariables(mean_left, var_left, log_nf_left, left);


                // update mu, v and z for the right branch
                double mean_right = mu_right[i];
                double log_nf_right = p_right[i];
                double var_right = v_right[i];

                propagateAuxiliaryVariables(mean_right, var_right, log_nf_right, right);


                // merging rule
                // D_node(y) = D_left(y) * D_right(y)
                
                // mean
                double mean_node = (mean_left*var_right + mean_right*var_left) / (var_left+var_right);
                mu_node[i] = mean_node;

                // var
                double var_node = (var_left*var_right) / (var_left+var_right);
                v_node[i] = var_node;

                if ( use_missing_data == true )
                {
                    missing_node[i] = false;
                }
                
                // log_nf
                double contrast = mean_left - mean_right;
                double a = -1.0 * contrast * contrast / ( 2.0 *(var_left+var_right) );
                double b = 0.5 * log( 2*RbConstants::PI*(var_left+var_right) );
                double log_norm_factor = log_nf_left + log_nf_right + a - b;
                double lnl_node = log_norm_factor;
                p_node[i] = lnl_node; //+ p_left[i] + p_right[i];
                
            } // end-if we had missing states for subtrees
            
            if ( node.isRoot() == true )
            {
                // this pruning algorithm is 100% equivalent to the likelihood obtained
                // using generalized least squares with a variance-covariance matrix
                // for the residuals (r = y - theta), also called the vcv-method (introduced in Hansen 1997)
                double root_state = computeRootState();
                p_node[i] += RbStatistics::Normal::lnPdf( root_state, v_node[i], mu_node[i]);

            }
            
        } // end for-loop over all characters
        
        
    } // end if we need to compute something for this node.
    
}



void PhyloOrnsteinUhlenbeckPruning::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already invalid, then all ancestral nodes must have been invalidated too
    if ( node_likelihoods.is_valid(index) )
    {
        // the root doesn't have an ancestor
        if ( !n.isRoot() )
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
void PhyloOrnsteinUhlenbeckPruning::invalidateBranchAndAncestors( const TopologyNode &n )
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
 * Invalidate all recomputed pruning entries while leaving fixed tip observations clean.
 * This is the full-cache fallback for topology and global parameter changes.
 */
void PhyloOrnsteinUhlenbeckPruning::invalidateInternalNodes( void )
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


void PhyloOrnsteinUhlenbeckPruning::resetValue()
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
    
    
    use_missing_data = false;
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();

    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            size_t index = (*it)->getIndex();
            ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );

            NodeLikelihood &tip_cache = node_likelihoods.init_for_writing(index);
            tip_cache.partial_likelihoods.assign(this->num_sites, 0.0);
            tip_cache.means.resize(this->num_sites);
            tip_cache.variances.assign(this->num_sites, 0.0);
            tip_cache.missing_data.assign(this->num_sites, false);

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


void PhyloOrnsteinUhlenbeckPruning::restoreSpecialization(void)
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.restore();
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
    size_t root_index = root.getIndex();
    
    // get the pointers to the partial likelihoods of the left and right subtree
    const std::vector<double> &p_root = this->node_likelihoods[root_index].partial_likelihoods;
    
    // sum the log-likelihoods for all sites together
    double sum_partial_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        sum_partial_probs += p_root[site];
    }
    
    return sum_partial_probs;
}


void PhyloOrnsteinUhlenbeckPruning::snapshotSpecialization( void )
{
    node_likelihoods.snapshot();
}


/*
 * Mark OU pruning likelihood cache entries dirty after a dependency changes.
 * Snapshot state is stored by IndexedSnapshotCache before this invalidation hook runs.
 */
void PhyloOrnsteinUhlenbeckPruning::invalidateSpecialization( const DagNode* affecter, bool fullyInvalidateSelf )
{
    // Branch-parameter changes can invalidate targeted branches; other changes fall back to full local invalidation.
    const TypedDagNode< RbVector< double > > *branch_parameter = NULL;
    if ( affecter == this->heterogeneous_alpha )
    {
        branch_parameter = this->heterogeneous_alpha;
    }
    else if ( affecter == this->heterogeneous_sigma )
    {
        branch_parameter = this->heterogeneous_sigma;
    }
    else if ( affecter == this->heterogeneous_theta )
    {
        branch_parameter = this->heterogeneous_theta;
    }

    if ( branch_parameter != NULL )
    {
        const std::set<size_t> &indices = branch_parameter->getTouchedElementIndices();
        
        // maybe all elements changed or the touched-element flags were not set precisely
        if ( indices.size() == 0 )
        {
            // just flag everything for recomputation
            fullyInvalidateSelf = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                this->invalidateBranchAndAncestors( *nodes[*it] );
            }
        }
    }

    if ( affecter == this->root_state )
    {
        recursivelyFlagNodeDirty( this->tau->getValue().getRoot() );
    }

    if ( affecter == this->homogeneous_alpha || affecter == this->homogeneous_sigma || affecter == this->homogeneous_theta )
    {
        fullyInvalidateSelf = true;
    }
    else if ( branch_parameter == NULL && affecter != this->root_state && affecter != this->tau )
    {
        fullyInvalidateSelf = true;
    }

    if ( affecter == this->dag_node )
    {
        resetValue();
    }
    
    if ( fullyInvalidateSelf )
    {
        invalidateInternalNodes();
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
