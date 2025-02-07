#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

#include "DistributionDirichlet.h"
#include "EmpiricalTreeTopologyProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TreeUtilities.h"
#include "Cloner.h"
#include "Proposal.h"
#include "RbOrderedSet.h"
#include "RbVector.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "TraceTree.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
EmpiricalTreeTopologyProposal::EmpiricalTreeTopologyProposal( StochasticNode<Tree> *t, const TraceTree &tt, double a ) : Proposal(),
    tree( t ),
    trees( tt.getValues() ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void EmpiricalTreeTopologyProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
EmpiricalTreeTopologyProposal* EmpiricalTreeTopologyProposal::clone( void ) const
{
    
//    return new EmpiricalTreeTopologyProposal( *this );
    return NULL;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& EmpiricalTreeTopologyProposal::getProposalName( void ) const
{
    static std::string name = "EmpiricalTopology";
    
    return name;
}


/**
 * Perform the proposal.
 *
 *
 *
 * \return The hastings ratio.
 */
double EmpiricalTreeTopologyProposal::doProposal( void )
{
    
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    double tree_age = tau.getRoot().getAge();
    
    size_t index = size_t( rng->uniform01() * trees.size() ) ;
    
    Tree *proposed_tree = trees[index].clone();
    TopologyNode *new_root = new TopologyNode( proposed_tree->getNumberOfNodes()+1 );
    
    TopologyNode &old_root = proposed_tree->getRoot();
    TopologyNode &og = old_root.getChild(2);
    old_root.removeChild(&og);
    new_root->addChild(&old_root);
    new_root->addChild( &og );
    old_root.setParent( new_root );
    og.setParent( new_root );
    
    double midpoint = proposed_tree->getBranchLengthForNode(og) / 2.0;
    old_root.setBranchLength( midpoint );
    og.setBranchLength( midpoint );
    
    proposed_tree->setRoot( new_root, true);
    proposed_tree->setRooted( true );
    
    // now we need to randomly root the tree
    const std::vector<TopologyNode*> nodes = proposed_tree->getNodes();
    size_t n_tips = proposed_tree->getNumberOfTips();
    size_t node_index = n_tips + size_t(floor(rng->uniform01() * (nodes.size()-n_tips)) );
    
    TopologyNode &random_root = *nodes[node_index];
    proposed_tree->reroot(random_root, false, true);
    
    std::vector<double> a = std::vector<double>(nodes.size()-n_tips-1,alpha);
    std::vector<double> proposed_ages = RbStatistics::Dirichlet::rv(a, *rng);
    
    stored_root_age = tree_age;
    RbOrderedSet<const TopologyNode*> current_lineages;
    RbOrderedSet<const TopologyNode*> available_parents;
    for (size_t i=n_tips; i<nodes.size(); ++i)
    {
        const TopologyNode &n = tau.getNode(i);
        current_lineages.insert( &n );
        
        if ( n.isRoot() == false )
        {
            const TopologyNode &parent = n.getParent();
            const TopologyNode *sibling = &parent.getChild(0);
            if ( sibling == &n )
            {
                sibling = &parent.getChild(1);
            }
            
            if ( current_lineages.find(sibling) != current_lineages.end() )
            {
                available_parents.insert( &parent );
            }
            
        }
        stored_ages[i-n_tips] = n.getAge() / tree_age;
    }

    double forward_prob = RbStatistics::Dirichlet::lnPdf(a, proposed_ages);
    double backward_prob = RbStatistics::Dirichlet::lnPdf(a, stored_ages);
    
    while ( available_parents.size() > 1 )
    {
        backward_prob += log( available_parents.size() );
        
        // I need to find the minimum age
        
    }
    
    std::sort( proposed_ages.begin(), proposed_ages.end() );
    std::set<TopologyNode*> tips;
    for (size_t i=0; i<n_tips; ++i)
    {
        TopologyNode &n = proposed_tree->getTipNode(i).getParent();
        
        bool viable = true;
        viable = viable && (n.isRoot() == false);
        for (size_t j=0; j<n.getNumberOfChildren(); ++j)
        {
            viable = viable && n.getChild(j).isTip();
        }
        if ( viable == true )
        {
            tips.insert( &n );
        }
        
    }
    size_t age_index = 0;
    while ( age_index < proposed_ages.size() )
    {
        size_t tip_index = size_t(floor(rng->uniform01() * tips.size()) );
        std::set<TopologyNode*>::const_iterator it(tips.begin());
        
        // 'advance' the iterator x times
        advance(it,tip_index);
        
    }
    
    // compute the Hastings ratio
    double lnHastingsratio = backward_prob - forward_prob;
    
    return lnHastingsratio;
}


/**
 *
 */
void EmpiricalTreeTopologyProposal::prepareProposal( void )
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
void EmpiricalTreeTopologyProposal::printParameterSummary(std::ostream &o) const
{
    
    o << "alpha = " << alpha;
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void EmpiricalTreeTopologyProposal::undoProposal( void )
{
    
    Tree& tau = tree->getValue();
    
    TopologyNode& node = tau.getRoot();
    
    // undo the proposal
    TreeUtilities::setAges(node, stored_ages );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void EmpiricalTreeTopologyProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
    
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void EmpiricalTreeTopologyProposal::tune( double rate )
{
    
    if ( rate > 0.234 )
    {
        alpha *= (1.0 + ((rate-0.234)/0.766) );
    }
    else
    {
        alpha /= (2.0 - rate/0.234 );
    }
    
}

