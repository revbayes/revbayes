#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "DistributionNormal.h"
#include "PhyloBrownianCharacterHistoryProcess.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "AbstractPhyloContinuousCharacterHistoryProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloBrownianCharacterHistoryProcess::PhyloBrownianCharacterHistoryProcess(const TypedDagNode<Tree> *t, size_t ns) :
    AbstractPhyloContinuousCharacterHistoryProcess( t, ns ),
    partial_likelihoods( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) )
{
    
    // We don't want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // now we need to reset the value
    this->redrawValue();
    
    // we need to reset the internal states
    resetValue();
    
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloBrownianCharacterHistoryProcess::~PhyloBrownianCharacterHistoryProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloBrownianCharacterHistoryProcess* PhyloBrownianCharacterHistoryProcess::clone( void ) const
{
    return new PhyloBrownianCharacterHistoryProcess( *this );
}


double PhyloBrownianCharacterHistoryProcess::computeLnProbability( void )
{
    
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        dirty_nodes = std::vector<bool>(tau->getValue().getNumberOfNodes(), true);
    }
    
//    // compute the ln probability by recursively calling the probability calculation for each node
//    const TopologyNode &root = this->tau->getValue().getRoot();
//    
//    // we start with the root and then traverse down the tree
//    size_t rootIndex = root.getIndex();
//    
//    // only necessary if the root is actually dirty
//    if ( this->dirty_nodes[rootIndex] )
//    {
//        
//        recursiveComputeLnProbability( root, rootIndex );
//        
//        // start by filling the likelihood vector for the children of the root
//        if ( root.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
//        {
////            recursiveComputeLnProbability( root, rootIndex );
//        }
//        else if ( root.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
//        {
////            const TopologyNode &left = root.getChild(0);
////            size_t left_index = left.getIndex();
////            recursiveComputeLnProbability( left, left_index );
////            const TopologyNode &right = root.getChild(1);
////            size_t right_index = right.getIndex();
////            recursiveComputeLnProbability( right, right_index );
////            const TopologyNode &middle = root.getChild(2);
////            size_t middleIndex = middle.getIndex();
////            recursiveComputeLnProbability( middle, middleIndex );
//            
//            //            computeRootLikelihood( rootIndex, left_index, right_index, middleIndex );
//            
//        }
//        else
//        {
//            throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
//        }
//        
//        // sum the partials up
//        this->lnProb = sumRootLikelihood();
//        
//    }
    
    return this->ln_prob;
    
}

void PhyloBrownianCharacterHistoryProcess::computeNodeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    if ( node.isRoot() == false && dirty_nodes[node_index] )
    {
        
        // get the branch length
        
        // get the branch rate
        
        for (int i=0; i < this->num_sites; ++i)
        {
            // get the site rate
            
            // get the ancestral state
            
            // get the current state

            // compute the transition probability
        }
        
        
        
    }
    
}


void PhyloBrownianCharacterHistoryProcess::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    flagNodeDirty( n );
    
}

void PhyloBrownianCharacterHistoryProcess::flagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirty_nodes[index] )
    {
        
        // set the flag
        dirty_nodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_nodes[index] == false )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }
        
    }
    
    // TODO: also need to mark the descendant nodes as dirty
    
}

void PhyloBrownianCharacterHistoryProcess::keepSpecialization( const DagNode* affecter )
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


//void PhyloBrownianCharacterHistoryProcess::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
//{
//
//    // check for recomputation
//    if ( node.isTip() == false && dirty_nodes[node_index] )
//    {
//        // mark as computed
//        dirty_nodes[node_index] = false;
//
//        std::vector<double> &p_node  = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
//        std::vector<double> &mu_node  = this->contrasts[this->active_likelihood[node_index]][node_index];
//
//        
//        // get the number of children
//        size_t num_children = node.getNumberOfChildren();
//        
//        for (size_t j = 1; j < num_children; ++j)
//        {
//        
//            size_t left_index = node_index;
//            const TopologyNode *left = &node;
//            if ( j == 1 )
//            {
//                left = &node.getChild(0);
//                left_index = left->getIndex();
//                recursiveComputeLnProbability( *left, left_index );
//            }
//            
//            const TopologyNode &right = node.getChild(j);
//            size_t right_index = right.getIndex();
//            recursiveComputeLnProbability( right, right_index );
//
//            const std::vector<double> &p_left  = this->partial_likelihoods[this->active_likelihood[left_index]][left_index];
//            const std::vector<double> &p_right = this->partial_likelihoods[this->active_likelihood[right_index]][right_index];
//
//            // get the per node and site contrasts
//            const std::vector<double> &mu_left  = this->contrasts[this->active_likelihood[left_index]][left_index];
//            const std::vector<double> &mu_right = this->contrasts[this->active_likelihood[right_index]][right_index];
//
//            // get the propagated uncertainties
//            double delta_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
//            double delta_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];
//
//            // get the scaled branch lengths
//            double v_left  = 0;
//            if ( j == 1 )
//            {
//                v_left = this->computeBranchTime(left_index, left->getBranchLength());
//            }
//            double v_right = this->computeBranchTime(right_index, right.getBranchLength());
//
//            // add the propagated uncertainty to the branch lengths
//            double t_left  = v_left  + delta_left;
//            double t_right = v_right + delta_right;
//
//            // set delta_node = (t_l*t_r)/(t_l+t_r);
//            this->contrast_uncertainty[this->active_likelihood[node_index]][node_index] = (t_left*t_right) / (t_left+t_right);
//
//            double stdev = sqrt(t_left+t_right);
//            for (int i=0; i<this->num_sites; i++)
//            {
//
//                mu_node[i] = (mu_left[i]*t_right + mu_right[i]*t_left) / (t_left+t_right);
//
//                // get the site specific rate of evolution
//                double standDev = this->computeSiteRate(i) * stdev;
//
//                // compute the contrasts for this site and node
//                double contrast = mu_left[i] - mu_right[i];
//
//                // compute the probability for the contrasts at this node
//                double lnl_node = RbStatistics::Normal::lnPdf(0, standDev, contrast);
//
//                // sum up the probabilities of the contrasts
//                p_node[i] = lnl_node + p_left[i] + p_right[i];
//
//            } // end for-loop over all sites
//
//        } // end for-loop over all children
//        
//    } // end if we need to compute something for this node.
//
//}

void PhyloBrownianCharacterHistoryProcess::resetValue( void )
{
    
    // make sure we have internal node states
    simulateInternalNodeStates();
    
    // check if the vectors need to be resized
    partial_likelihoods = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );

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


void PhyloBrownianCharacterHistoryProcess::restoreSpecialization( const DagNode* affecter )
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

void PhyloBrownianCharacterHistoryProcess::simulateInternalNodeStates( void )
{
 
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxa = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    // simulate the root sequence
    ContinuousTaxonData &root = taxa[ tau->getValue().getRoot().getIndex() ];
    
    std::vector<double> root_states = simulateRootCharacters(num_sites);
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // create the character
        double c = root_states[i];
        
        // add the character to the sequence
        root.addCharacter( c );
    }
    
    this->value->addTaxonData( root );
    
    // recursively simulate the sequences
    simulateRecursivelyInternal( tau->getValue().getRoot(), taxa );
    
}


void PhyloBrownianCharacterHistoryProcess::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        TopologyNode &child = *(*it);
        
        // name the taxon if it isn't already
        if ( child.getTaxon().getName() == "" )
        {
            child.setName("node " + StringUtilities::to_string( child.getIndex() ) );
        }
        
        // get the branch length for this child
        double branch_length = tau->getValue().getBranchLengthForNode(child);
        
        // get the branch specific rate
        double branch_sigma = sqrt( computeBranchTime( child.getIndex(), branch_length ) );
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double standDev = branch_sigma * computeSiteRate( i );
            
            // create the character
            double c = RbStatistics::Normal::rv( parent_state, standDev, *rng);
            
            // add the character to the sequence
            taxon.addCharacter( c );
        }
        
        taxon.setTaxon( child.getTaxon() );
        if ( child.isTip() == false )
        {
            simulateRecursively( child, taxa );
        }
        
    }
    
}

void PhyloBrownianCharacterHistoryProcess::simulateRecursivelyInternal( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        
        if ( child.isTip() == false )
        {

            // get the branch length for this child
            double branch_length = tau->getValue().getBranchLengthForNode(child);
            
            // get the branch specific rate
            double branch_sigma = sqrt( computeBranchTime( child.getIndex(), branch_length ) );
            
            ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
            for ( size_t i = 0; i < num_sites; ++i )
            {
                // get the ancestral character for this site
                double parent_state = parent.getCharacter( i );
                
                // compute the standard deviation for this site
                double standDev = branch_sigma * computeSiteRate( i );
                
                // create the character
                double c = RbStatistics::Normal::rv( parent_state, standDev, *rng);
                
                // add the character to the sequence
                taxon.addCharacter( c );
            }
            
            taxon.setTaxon( child.getTaxon() );
            
            // add the taxon data to the character data
            this->value->addTaxonData( taxon );
            
            // recursively simulate the sequences
            simulateRecursively( child, taxa );
            
        }
        
    }
    
}

void PhyloBrownianCharacterHistoryProcess::touchSpecialization( const DagNode* affecter, bool touchAll )
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
                this->flagNodeDirty( *nodes[*it] );
            }
        }
    }
    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
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
void PhyloBrownianCharacterHistoryProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloContinuousCharacterHistoryProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else
    {
        AbstractPhyloContinuousCharacterHistoryProcess::swapParameterInternal(oldP, newP);
    }
    
}

