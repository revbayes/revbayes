#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "PhyloOrnsteinUhlenbeckThreePoint.h"
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

PhyloOrnsteinUhlenbeckThreePoint::PhyloOrnsteinUhlenbeckThreePoint(const TypedDagNode<Tree> *t, size_t ns) : AbstractPhyloContinuousCharacterProcess( t, ns ),
partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
contrasts( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
contrast_uncertainty( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
normalizing_constants( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 1.0) ) ) ),
active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
changed_nodes( std::vector<bool>(this->num_nodes, false) ),
dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
obs()
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
    
    // we need to reset the contrasts
    resetValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloOrnsteinUhlenbeckThreePoint::~PhyloOrnsteinUhlenbeckThreePoint( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloOrnsteinUhlenbeckThreePoint* PhyloOrnsteinUhlenbeckThreePoint::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckThreePoint( *this );
}


double PhyloOrnsteinUhlenbeckThreePoint::computeBranchAlpha(size_t branch_idx) const
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


double PhyloOrnsteinUhlenbeckThreePoint::computeBranchSigma(size_t branch_idx) const
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


double PhyloOrnsteinUhlenbeckThreePoint::computeBranchTheta(size_t branch_idx) const
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


double PhyloOrnsteinUhlenbeckThreePoint::computeRootState( void ) const
{
    
    // get the root-state parameter
    double root_state = this->root_state->getValue();
    
    return root_state;
}


double PhyloOrnsteinUhlenbeckThreePoint::computeLnProbability( void )
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
    
    size_t num_tips = tau->getValue().getNumberOfTips();
    
    double alpha    = homogeneous_alpha->getValue();
    double theta    = homogeneous_theta->getValue();
    double rs       = root_state->getValue();
    double sigma    = homogeneous_sigma->getValue();
    double sigma2   = sigma*sigma;
    
    // we need to root to tip distance, which is simply the root age for ultrametric trees
    double root_age = root.getAge();
    
    std::vector<double> X = std::vector<double>(num_tips, 0.0);
    for (size_t i=0; i<num_tips; ++i)
    {
        // need to change this for non-ultrametric trees
        // X[i] = rs *exp(alpha*dis[i]) + theta*(1-exp(alpha*dis[i]));
        X[i] = rs *exp(alpha*root_age) + theta*(1-exp(alpha*root_age));
    }
    
    
//    tree = transf.branch.lengths(phy,model,parameters=prm,
//                                 check.pruningwise=F,check.ultrametric=F,D=D,check.names=F)$tree
    
    std::vector<double >branches_lengths = transformBranchLengths();
    
//    tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
//           as.integer(1),as.integer(d),as.integer(ROOT),as.double(tree$root.edge),as.double(tree$edge.length),
//           as.integer(des), as.integer(anc), as.double(as.vector(y)), as.double(as.vector(X)),
    std::vector<double> output;
    threePoint(output, branches_lengths, obs[0], X);

//## logdetV, 1'V^{-1}1, y'V^{-1}1, y'V^{-1}y, X'V^{-1}1, X'V^{-1}X, X'V^{-1}y
//    comp = list(vec11=tmp[2], y1=tmp[3], yy=tmp[4], X1=tmp[5:(4+d)],
//                XX=matrix(tmp[(5+d):(ole-d)], d,d),Xy=tmp[(ole-d+1):ole],logd=tmp[1])
    size_t d = 1;
    size_t ole = 4 + 2*d + d*d; // # output length
    double logd = output[0];
    double yy   = output[3];
    double xx   = output[4+d];
    double xy   = output[(ole-d)];
    double lnl = num_tips*log(RbConstants::TwoPI);
    lnl += logd;
    lnl += num_tips*log(sigma2/2.0/alpha);
    lnl += 2.0*alpha/sigma2*(yy - 2.0*xy + xx);
    lnl /= -2.0;
    
    return lnl;
}



void PhyloOrnsteinUhlenbeckThreePoint::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


void PhyloOrnsteinUhlenbeckThreePoint::keepSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckThreePoint::recursivelyFlagNodeDirty( const TopologyNode &n )
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


void PhyloOrnsteinUhlenbeckThreePoint::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrasts = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrast_uncertainty = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );
    normalizing_constants = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 1.0) ) );
    
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
    
    size_t num_species = this->tau->getValue().getNumberOfTips();
    obs = std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_species, 0.0) );
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(site_indices[site]);
                obs[site][(*it)->getIndex()] = c;
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


void PhyloOrnsteinUhlenbeckThreePoint::restoreSpecialization( const DagNode* affecter )
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


void PhyloOrnsteinUhlenbeckThreePoint::setAlpha(const TypedDagNode<double> *a)
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


void PhyloOrnsteinUhlenbeckThreePoint::setAlpha(const TypedDagNode<RbVector<double> > *a)
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


void PhyloOrnsteinUhlenbeckThreePoint::setRootState(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckThreePoint::setSigma(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckThreePoint::setSigma(const TypedDagNode<RbVector<double> > *s)
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


void PhyloOrnsteinUhlenbeckThreePoint::setTheta(const TypedDagNode<double> *t)
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


void PhyloOrnsteinUhlenbeckThreePoint::setTheta(const TypedDagNode<RbVector<double> > *t)
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


void PhyloOrnsteinUhlenbeckThreePoint::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        double branch_length = tau->getValue().getBranchLengthForNode(child);
        
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


std::vector<double> PhyloOrnsteinUhlenbeckThreePoint::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState();
    }
    
    return chars;
}


double PhyloOrnsteinUhlenbeckThreePoint::sumRootLikelihood( void )
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


void PhyloOrnsteinUhlenbeckThreePoint::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_sigma )
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


std::vector<double> PhyloOrnsteinUhlenbeckThreePoint::transformBranchLengths( void ) const
{
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root                    = this->tau->getValue().getRoot();
    const std::vector<TopologyNode*> &nodes     = this->tau->getValue().getNodes();
    
    size_t num_tips     = tau->getValue().getNumberOfTips();
    size_t num_internal = tau->getValue().getNumberOfInteriorNodes()+1;
    size_t num_branches = this->tau->getValue().getNumberOfNodes();
    
    double alpha    = homogeneous_alpha->getValue();
//    double theta    = homogeneous_theta->getValue();
//    double rs       = root_state->getValue();
//    double sigma    = homogeneous_sigma->getValue();
//    double sigma2   = sigma*sigma;
    
    double root_age = root.getAge();
    
    std::vector<double> times = std::vector<double>(num_internal, 0.0);
    for (size_t i=0; i<nodes.size(); ++i)
    {
        TopologyNode *the_node = nodes[i];
        if ( the_node->isTip() == false )
        {
            size_t index = the_node->getIndex() - num_tips;
            if ( index >= times.size() )
                throw RbException("out of bounds");
            times[ index ] = the_node->getAge();
        }
    }
//    errEdge = errEdge*exp(-2*alpha*D[des[externalEdge]]) # adjust measurement errors for OU models

    std::vector<double> dist_from_root = std::vector<double>(num_internal, 0.0);
    for (size_t i=0; i<num_internal; ++i)
    {
        if ( i >= dist_from_root.size() )
            throw RbException("out of bounds");
        dist_from_root[i] = exp(-2*alpha*times[i]) * (1 - exp(-2*alpha*(root_age-times[i])));
    }
    
    std::vector<double> anc_times   = std::vector<double>(num_branches, 0.0);
    for (size_t i=0; i<nodes.size(); ++i)
    {
        TopologyNode *the_node = nodes[i];
        if ( the_node->isRoot() == false )
        {
            if ( the_node->getIndex() >= anc_times.size() )
                throw RbException("out of bounds");
            anc_times[ the_node->getIndex() ] = dist_from_root[ the_node->getParent().getIndex() - num_tips ];
        }
    }
    
    std::vector<double> d2 = std::vector<double>(num_branches, 0.0);
    for (size_t i=0; i<num_branches; ++i)
    {
        TopologyNode *the_node = nodes[i];
        size_t node_index = the_node->getIndex();
        if ( the_node->isRoot() == false && the_node->isTip() == true )
        {
            if ( node_index >= d2.size() )
                throw RbException("out of bounds");
            d2[node_index] = (1-exp(-2*alpha*root_age));
            // need to adapt this for non-ultramtric trees
            // d2[node_index] = exp(-2*alpha*D[des[externalEdge]]) * (1-exp(-2*alpha*(Tmax-D[des[externalEdge]])))
        }
        else if ( the_node->isRoot() == false && the_node->isInternal() == true )
        {
            if ( (node_index-num_tips) >= dist_from_root.size() || node_index >= d2.size() )
                throw RbException("out of bounds");
            d2[node_index] = dist_from_root[node_index-num_tips];
//            d2[node_index] = dist_from_root[des[!externalEdge]-n]
        }
        
    }
    
    std::vector<double> edge_lengths = std::vector<double>(num_branches, 0.0);
    for (size_t i=0; i<num_branches; ++i)
    {
        if ( i >= edge_lengths.size() || i >= d2.size() || i >= anc_times.size() )
            throw RbException("out of bounds");
        edge_lengths[i] = d2[i] - anc_times[i];
    }
    
    return edge_lengths;
//        root.edge = min(distFromRoot)
//        diagWeight = exp(alpha*D)
//
//        edge.length[externalEdge] = edge.length[externalEdge] + errEdge # add measurement errors to the tree
//        phy$edge.length = edge.length
//        phy$root.edge = root.edge

}


/** Swap a parameter of the distribution */
void PhyloOrnsteinUhlenbeckThreePoint::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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


void PhyloOrnsteinUhlenbeckThreePoint::threePoint(std::vector<double> &output, const std::vector<double> &bl, const std::vector<double> &obs, const std::vector<double> &mean)
{
    
    const std::vector<TopologyNode*> &nodes = this->tau->getValue().getNodes();
    
    
    int n   = this->tau->getValue().getNumberOfTips();              // number of tips
    int N   = 2*n-2;                                                // number of edges, excluding root edge
    int pN  = this->tau->getValue().getNumberOfInteriorNodes()+1;   // number of internal nodes
    int npN = n + pN;                                               // total # nodes, including root node. N+1 if binary tree.
    int dY  = 1;                                                    // number of columns in Y (dependent variables, usually dY=1)
    int dX  = 1;                                                    // number of columns in X (independent variables)
    int r   = this->tau->getValue().getRoot().getIndex();           // root node ID (in R)
    
    double rootEdge= 0; // length of root edge
    // transb: vector of edge lengths
    // des:    vector of descendents (one for each edge)
    // anc:    vector of ancestor nodes
    // y:  vector of response values
    // X:  vectorized matrix X, where matrix model is: y = X*beta + error
    // output: vector of size 4+2dX+dX^2: logdet(V), 1'V^{-1}1,
    //         y'V^{-1}1 (dY values), y'V^{-1}y (dY^2 values),
    //         X'V^{-1}1 (dX values), X'V^{-1}X (dX^2 values), X'V^{-1}y (dX*dY values)
    
    std::vector<double> logd    = std::vector<double>(npN,      0.0);
    std::vector<double> yy      = std::vector<double>(npN*dY*dY,0.0);
    std::vector<double> xy      = std::vector<double>(npN*dX*dY,0.0);
    std::vector<double> xx      = std::vector<double>(npN*dX*dX,0.0);
    
    // logd[i-1] = log|V_subtree| for subtree at *node* i, or parent edge of node i.
    std::vector<double> vec11= std::vector<double>(npN,0.0);
    // first (Y'V^{-1}Y)[1,1] at all nodes, then (Y'V^{-1}Y)[2,1] at all nodes, etc.
    std::vector<double> y1 = std::vector<double>(npN*dY,0.0);
    std::vector<double> X1 = std::vector<double>(npN*dX,0.0);
    std::vector<int> zero =std::vector<int>(npN,0);
    // zero[i-1] will be -1 if node i has no children edge of length 0
    //                    d if node i has exactly 1 child edge of length 0,
    //                      to child node d+1.
    // error returned if some node has 2 or more child edges of length 0
    for (int iedge=0; iedge<N+1; iedge++)
        zero[iedge] = -1;
    
    // loop over all N+1 edges, with root edge last. Assumes postorder traversal.
    for (int iedge=0; iedge < N+1; iedge++)
    {
        double el, invel;  // edge length and its inverse
        int di, anci=0;
        if ( iedge >= nodes.size() )
            throw RbException("out of bounds");
        TopologyNode *the_node = nodes[iedge];
        
        if (iedge<N){      // non-root edge
            el      = bl[iedge];
            di      = the_node->getIndex();// subtree at node di+1.
            anci    = the_node->getParent().getIndex();
        } else {           // iedge=N -> root edge
            el = rootEdge;
            di = r; // but anci meaningless
        }
        
        int iY1 = di; // index in y1 for current edge and column j in Y: di + npN*j
        int iX1 = di; //          X1                                  X: di + npN*j
        int iYY = di; // index in yy for columns k and j in Y: di + npN*k + npN*dY*j
        int iXY = di; // index in Xy for column k in X and j in Y: di + npN*k + npN*dX*j
        int iXX = di; // di + npN*k + npN*dX*j for XX: columns k and j in X.
        
        if (di<n){ // external edge
            if(el>0) invel = 1/el;// y1, X1 will actually hold ybar*p and Xbar*p.
            else     invel = 1.0; // in this case, we need the staight y'1, y'y, X'1, X'X etc.
            if (el>0){
                logd[di]=log(el);
                vec11[di]=invel;
            } else{
                if (zero[anci] >= 0) // anci already found with 1 child edge of length 0.
                    throw RbException("two or more sister external edges have length 0, V is singular\n (node in pruning-wise-ordered tree)");
                else
                    zero[anci] = di; // which is >= 0
            }
            int jY  = di; // index in y for observations in species di, column j: di + n*j
            for (int j=0; j<dY; j++){
                y1[iY1] = obs[jY]*invel;
                int kY = di; // index in y for observation in species di, column k: di + n*k
                for (int k=0; k<dY; k++){
                    yy[iYY] = y1[iY1] * obs[kY];
                    iYY += npN;
                    kY  += n;
                }
                int kX = di; // index in X for observation in species di, column k: di + n*k
                for (int k=0; k<dX; k++){
                    xy[iXY] = y1[iY1] * mean[kX];
                    iXY += npN;
                    kX  += n;
                }
                iY1 += npN;
                jY  += n;
            }
            int jX  = di; // di + n*j   for X, species di, column j.
            for (int j=0; j<dX; j++){
                X1[iX1] = mean[jX] * invel;
                int kX = di;
                for (int k=0; k<dX; k++){
                    xx[iXX] = X1[iX1] * mean[kX];
                    iXX += npN;
                    kX  += n;
                }
                iX1 += npN;
                jX  += n;
            }
        }
        else{ // internal edge. contributions from descendants of di have already been collected.
            int goodchildren = 1;
            double ev, ev2;
            if (zero[di] >= 0) { // exactly 1 child edge of length 0, to descendant d0=zero[di]
                if (el<=0)
                    throw RbException("One external edge and its parent both have length 0\n (node in pruning-wise-ordered tree). To avoid this situation,\n please make a polytomy or resolve it differently ");
                goodchildren = 0;
                // still assuming di has more than 1 child, with current vec11[di]>0.
            }
            if (goodchildren) {
                logd[di] += log(1+ el * vec11[di]);
                ev  = el/(1+ el * vec11[di]);
                ev2 =  1/(1+ el * vec11[di]);
                for (int j=0; j<dY; j++){
                    double tmp1 = ev * y1[iY1];
                    int kY1 = di; // index in ybar, column k: di + npN*k
                    for (int k=0; k<dY; k++){
                        yy[iYY] -= tmp1 * y1[kY1];
                        iYY += npN;
                        kY1 += npN;
                    }
                    int kX1 = di; // index in Xbar, column k: di + npN*k
                    for (int k=0; k<dX; k++){
                        xy[iXY] -= tmp1 * X1[kX1];
                        iXY += npN;
                        kX1 += npN;
                    }
                    iY1 += npN;
                }
                for (int j=0; j<dX; j++){
                    double tmp1 = ev * X1[iX1];
                    int kX1 = di; // index in Xbar, column k: di + npN*k
                    for (int k=0; k<dX; k++){
                        xx[iXX] -= tmp1 * X1[kX1];
                        iXX += npN;
                        kX1 += npN;
                    }
                    iX1 += npN;
                }
            }
            else { // 1 bad child
                logd[di] += log(el);
                int d0 = zero[di];
                double fac = 1/el + vec11[di];
                // X'V^{-1}y = Xytilde + (1/t +ptilde_A) X'd0 yd0
                //            - ptilde_A X'd0 y1tilde - ptilde_A X1tilde' yd0
                // tilde: values collected from "normal" children
                // d0:    values from the special child with edge length 0
                iY1 = 0; iX1 = 0;
                for (int j=0; j<dY; j++){
                    double tmp1 = fac * y1[d0+iY1];
                    int kY1 = 0; // index in ybar, column k: (di or d0) + npN*k
                    for (int k=0; k<dY; k++){
                        yy[iYY] += tmp1 * y1[d0+kY1] - y1[d0+iY1]*y1[di+kY1] - y1[di+iY1]*y1[d0+kY1];
                        iYY += npN;
                        kY1 += npN;
                    }
                    int kX1 = 0; // index in Xbar, column k: (di or d0) + npN*k
                    for (int k=0; k<dX; k++){
                        xy[iXY] += tmp1 * X1[d0+kX1] - y1[d0+iY1]*X1[di+kX1] - y1[di+iY1]*X1[d0+kX1];
                        iXY += npN;
                        kX1 += npN;
                    }
                    iY1 += npN;
                }
                for (int j=0; j<dX; j++){
                    double tmp1 = fac * X1[d0+iX1];
                    int kX1 = 0;
                    for (int k=0; k<dX; k++){
                        xx[iXX] += tmp1 * X1[d0+kX1] - X1[d0+iX1]*X1[di+kX1] - X1[di+iX1]*X1[d0+kX1];
                        iXX += npN;
                        kX1 += npN;
                    }
                    iX1 += npN;
                }
            }
            if (goodchildren) {
                // rescaling to get vec11 correct ("p" instead of "p_A"),
                // and Y1 = ybar*p, X1 = bar*p instead of ybar*p_A etc.
                iY1 = di;
                for (int j=0; j<dY; j++){
                    y1[iY1] *= ev2;
                    iY1 += npN;
                }
                iX1 = di;
                for (int j=0; j<dX; j++){
                    X1[iX1] *= ev2;
                    iX1 += npN;
                }
                vec11[di] *= ev2;
            }
            else {
                invel = 1/el;
                iY1 = 0;
                for (int j=0; j<dY; j++){
                    y1[di+iY1] = y1[zero[di]+iY1]*invel;
                    iY1 += npN;
                }
                iX1 = 0;
                for (int j=0; j<dX; j++){
                    X1[di+iX1] = X1[zero[di]+iX1]*invel;
                    iX1 += npN;
                }
                vec11[di] = invel;
            }
        }
        // next: collect contribution of iedge to its parent, i.e. to ancestor node anci.
        // *except* if root edge (anci meaningless) or if external edge of zero length
        if ((iedge < N) && ((di>=n) || (el>0))){
            logd[anci] += logd[di];
            iY1 = 0; iX1 = 0;
            iYY = 0; iXX = 0; iXY = 0;
            for (int j=0; j<dY; j++){
                y1[anci+iY1] += y1[di+iY1];
                for (int k=0; k<dY; k++){
                    yy[anci + iYY] += yy[di + iYY];
                    iYY += npN;
                }
                for (int k=0; k<dX; k++){
                    xy[anci + iXY] += xy[di + iXY];
                    iXY += npN;
                }
                iY1 += npN;
            }
            for (int j=0; j<dX; j++){
                X1[anci+iX1] += X1[di+iX1];
                for (int k=0; k<dX; k++){
                    xx[anci + iXX] += xx[di + iXX];
                    iXX += npN;
                }
                iX1 += npN;
            }
            vec11[anci] += vec11[di];
        }
    }
    // combine all results at root edge into output
    output.resize(2+dY+dY*dY+dX+dX*dX+dX*dY);
    output[0]=logd[r];
    output[1]=vec11[r];
    int p=2; // # parameters already logged
    int ikXY = r;
    for (int j=0; j<dY; j++){
        output[p+j] = y1[ikXY];
        ikXY += npN;
    }
    p += dY;
    ikXY = r;
    for (int j=0; j<dY*dY; j++){
        output[p+j] = yy[ikXY];
        ikXY += npN;
    }
    p += dY*dY;
    ikXY = r;
    for (int j=0; j<dX; j++){
        output[p+j] = X1[ikXY];
        ikXY += npN;
    }
    p += dX;
    ikXY = r;
    for (int j=0; j<dX*dX; j++){
        output[p+j] = xx[ikXY];
        ikXY += npN;
    }
    p += dX*dX;
    ikXY = r;
    for (int j=0; j<dX*dY; j++){
        output[p+j] = xy[ikXY];
        ikXY += npN;
    }
    
}


