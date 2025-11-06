#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "MixtureDistribution.h"
#include "NarrowExchangeRateMatrixProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "DagNode.h"
#include "MemberObject.h"
#include "Proposal.h"
#include "RateGenerator.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDistribution.h"

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NarrowExchangeRateMatrixProposal::NarrowExchangeRateMatrixProposal( StochasticNode<Tree> *t, const std::vector<StochasticNode<RateGenerator> *> &rm ) : Proposal(),
    tree( t ),
    rate_matrices( rm )
{
    // tell the base class to add the node
    addNode( tree );
    
    for (size_t i=0; i<rate_matrices.size(); ++i)
    {
        addNode( rate_matrices[i] );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NarrowExchangeRateMatrixProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NarrowExchangeRateMatrixProposal* NarrowExchangeRateMatrixProposal::clone( void ) const
{
    
    return new NarrowExchangeRateMatrixProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NarrowExchangeRateMatrixProposal::getProposalName( void ) const
{
    static std::string name = "NarrowExchangeRateMatrix";
    
    return name;
}


double NarrowExchangeRateMatrixProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}



/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new simplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NarrowExchangeRateMatrixProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root and neithor a direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    TopologyNode& parent = node->getParent();
    TopologyNode& grandparent = parent.getParent();
    TopologyNode* uncle = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == &parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    TopologyNode* brother = &parent.getChild( 0 );
    // check if we got the correct child
    if ( brother == node )
    {
        brother = &parent.getChild( 1 );
    }
    
    // we need to work with the times
    double parent_age   = parent.getAge();
    double uncles_age   = uncle->getAge();
    
    if ( uncles_age < parent_age )
    {
        failed = false;
        
        // now we store all necessary values
        stored_choosen_node     = node;
        stored_uncle            = uncle;
        
        // compute the backwards probabilities
        double backwards = 0.0;
        backwards += lnProposalProbabilityRateMatrix( *node,    stored_index_node,    false );
        backwards += lnProposalProbabilityRateMatrix( parent,   stored_index_parent,  false );
        backwards += lnProposalProbabilityRateMatrix( *brother, stored_index_brother, false );
        backwards += lnProposalProbabilityRateMatrix( *uncle,   stored_index_uncle,   false );
        
        // now exchange the two nodes
        grandparent.removeChild( uncle );
        parent.removeChild( node );
        grandparent.addChild( node );
        parent.addChild( uncle );
        node->setParent( &grandparent );
        uncle->setParent( &parent );
        
        // update the branches and compute the forward probabilities
        double forwards = 0.0;
        forwards += lnProposalProbabilityRateMatrix( *uncle,   proposed_index_uncle,   true );
        forwards += lnProposalProbabilityRateMatrix( parent,   proposed_index_parent,  true );
        forwards += lnProposalProbabilityRateMatrix( *brother, proposed_index_brother, true );
        forwards += lnProposalProbabilityRateMatrix( *node,    proposed_index_node,    true );
        
        return backwards - forwards;
    }
    else
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
}


double NarrowExchangeRateMatrixProposal::lnProposalProbabilityRateMatrix(const TopologyNode &node, size_t &old_category, bool update)
{
    size_t node_index = node.getIndex();
    StochasticNode<RateGenerator> *rm = rate_matrices[node_index];
    
    // potential affected nodes for likelihood computation
    RbOrderedSet<DagNode *> affected;
    rm->initiateGetAffectedNodes( affected );
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    MixtureDistribution<RateGenerator>& dist = static_cast<MixtureDistribution<RateGenerator> &>( rm->getDistribution() );
    
    // get the number of categories
    size_t n = dist.getNumberOfMixtureElements();
    
    // create a vector for the weights
    std::vector<double> weights = std::vector<double>(n,0);
    double sum_of_weights = 0.0;
    double max_weight = RbConstants::Double::neginf;
    
    // get the current index
    old_category = dist.getCurrentIndex();
    
    for (size_t i=0; i<n; ++i)
    {
        // set our new value
        dist.setCurrentIndex( i );
        
        // flag for likelihood recomputation
        rm->touch();
        
        // compute the likelihood of the new value
        double prior_ratio = rm->getLnProbability();
        double likelihood_ratio = 0.0;
        for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
        {
            likelihood_ratio += (*it)->getLnProbability();
        }
        weights[i] = prior_ratio + likelihood_ratio;
        
        if (max_weight < weights[i])
        {
            max_weight = weights[i];
        }
        
    }
    
    // normalize weights
    for (size_t i=0; i<n; ++i)
    {
        weights[i] = exp(weights[i] - max_weight);
        sum_of_weights += weights[i];
    }
    
    size_t new_category = old_category;
    if ( update == true )
    {
        double u = rng->uniform01() * sum_of_weights;
        new_category = 0;
        do
        {
            u -= weights[new_category];
            ++new_category;
        } while ( u > 0.0 );
        --new_category;
    }
    
    // set our new value
    dist.setCurrentIndex( new_category );
    
    double ln_p = log( weights[new_category] / sum_of_weights );
    
    return ln_p;
}



/**
 *
 */
void NarrowExchangeRateMatrixProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void NarrowExchangeRateMatrixProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NarrowExchangeRateMatrixProposal::undoProposal( void )
{
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        TopologyNode& parent        = stored_uncle->getParent();
        TopologyNode& grandparent   = stored_choosen_node->getParent();
        TopologyNode& node          = *stored_choosen_node;
        TopologyNode& uncle         = *stored_uncle;
        TopologyNode* brother       = &parent.getChild( 0 );
        // check if we got the correct child
        if ( brother == &uncle )
        {
            brother = &parent.getChild( 1 );
        }


        // now exchange the two nodes
        grandparent.removeChild( stored_choosen_node );
        parent.removeChild( stored_uncle );
        grandparent.addChild( stored_uncle );
        parent.addChild( stored_choosen_node );
        stored_uncle->setParent( &grandparent );
        stored_choosen_node->setParent( &parent );
        
        undoRateMatrixProposal( parent,   stored_index_parent );
        undoRateMatrixProposal( node,     stored_index_node );
        undoRateMatrixProposal( uncle,    stored_index_uncle );
        undoRateMatrixProposal( *brother, stored_index_brother );
        
    }
    
}


void NarrowExchangeRateMatrixProposal::undoRateMatrixProposal(const RevBayesCore::TopologyNode &node, size_t i)
{
    size_t node_index = node.getIndex();
    StochasticNode<RateGenerator> *rm = rate_matrices[node_index];
    
    MixtureDistribution<RateGenerator>& dist = static_cast<MixtureDistribution<RateGenerator> &>( rm->getDistribution() );
        
    // set our old value
    dist.setCurrentIndex( i );

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NarrowExchangeRateMatrixProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    
    for (size_t i=0; i<rate_matrices.size(); ++i)
    {
        if ( oldN == rate_matrices[i] )
        {
            rate_matrices[i] = static_cast<StochasticNode<RateGenerator>* >(newN) ;
        }
        
    }
    
    
}


void NarrowExchangeRateMatrixProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void NarrowExchangeRateMatrixProposal::tune( double rate )
{
    
    // nothing to tune
    
}

