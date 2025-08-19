#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "DistributionUniform.h"
#include "CollapseExpandFossilBranchProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TypedDagNode.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
CollapseExpandFossilBranchProposal::CollapseExpandFossilBranchProposal( StochasticNode<Tree> *n, TypedDagNode<double> *o ) : Proposal(),
    tau( n ),
    origin( o )
{
    // tell the base class to add the node
    addNode( tau );
    addNode( origin );
    if(! tau->getDistribution().allowsSA()) throw RbException("Setup includes a move for sampled ancestors but the corresponding tree distribution doesn't allow sampled ancestors");
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CollapseExpandFossilBranchProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CollapseExpandFossilBranchProposal* CollapseExpandFossilBranchProposal::clone( void ) const
{
    
    return new CollapseExpandFossilBranchProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CollapseExpandFossilBranchProposal::getProposalName( void ) const
{
    static std::string name = "CollapseExpandFossilBranch";
    
    return name;
}


double CollapseExpandFossilBranchProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double CollapseExpandFossilBranchProposal::doProposal( void )
{
    
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree &t = tau->getValue();
    
    std::vector<TopologyNode*> fossils;
    for (size_t i = 0; i < t.getNumberOfTips(); ++i)
    {
        TopologyNode* node = &t.getNode(i);

        if ( node->isFossil() == true )
        {
            fossils.push_back(node);
        }

    }

    if ( fossils.empty() )
    {
        failed = true;
        return 0;
    }

    // pick a random fossil node
    double u = rng->uniform01();
    storedNode = fossils[ size_t(u*fossils.size()) ];
    
    double hr = 0;
    if ( storedNode->isSampledAncestorTip() == true )
    {
        hr += expandBranch( *storedNode );
    }
    else
    {
        hr += collapseBranch( *storedNode );
    }
    
    return hr;
    
}



/* Move a fossil tip (brl > 0) to be ancestral (brl =0)
 
 __|__              __|__
|     |            |     |
     q|___     -->       |
      |   |    -->       |
      |   |p            q|___p
     r|                 r|
 
 1. Pick a fossil among those with brl > 0 (prob = 1/m)
 2. Set brl = 0
 */
double CollapseExpandFossilBranchProposal::collapseBranch(TopologyNode &n)
{
    
    // Get the parent and sibling of the chosen node
    TopologyNode &parent = n.getParent();
    TopologyNode *sibling = &parent.getChild( 0 );
    if ( sibling == &n )
    {
        sibling = &parent.getChild( 1 );
    }
    
    // determine lower and upper bound of backward move
    double min_age = n.getTaxon().getAgeRange().getMax();
    double max_age = parent.getAge();
    if ( parent.isRoot() )
    {
        if (origin != NULL)
        {
            max_age = origin->getValue();
        }
    }
    else
    {
        max_age = parent.getParent().getAge();
    }
    
    // test that the max age is larger than the min age
    if ( max_age <= min_age || n.getAge() < sibling->getAge() )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // store the old age of the parent
    storedAge = parent.getAge();
    
    // set the age of the parent node equal to my age
    parent.setAge( n.getAge() );
    n.setSampledAncestor( true );
    
    // compute the Jacobian term
    double lnJacobian = - log(max_age - min_age);
    
    return lnJacobian;
}



/* Move an ancestral fossil (brl = 0) to fossil tip (brl > 0)
 
    __|__              __|__
   |     |            |     |
         |       -->       q|___
         |       -->        |   |
        q|___p              |   |p
        r|                 r|
 
 1. Pich a fossil among those with brl = 0 (prob = 1/k)
 2. Propose brl from a uniform(0, ?) distribution
 */
double CollapseExpandFossilBranchProposal::expandBranch(TopologyNode &n)
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // Get the parent and sibling of the chosen node
    TopologyNode &parent = n.getParent();
    TopologyNode *sibling = &parent.getChild( 0 );
    if ( sibling == &n)
    {
        sibling = &parent.getChild( 1 );
    }
    
    // determine lower and upper bound of backward move
    double min_age = n.getTaxon().getAgeRange().getMax();
    double max_age = parent.getAge();
    if ( parent.isRoot() )
    {
        if (origin != NULL)
        {
            max_age = origin->getValue();
        }
    }
    else
    {
        max_age = parent.getParent().getAge();
    }
    
    // test that the max age is larger than the min age
    if ( max_age <= min_age )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    // store the old age of the parent
    storedAge = parent.getAge();
    
    // draw the new age for the parent node
    double new_age = (max_age-min_age) * rng->uniform01() + min_age;
    
    // set the age of the parent node equal to the new age
    n.setSampledAncestor( false );
    parent.setAge( new_age );
    
    // compute the Jacobian term
    double lnJacobian = log(max_age - min_age);
    
    return lnJacobian;
}



/**
 *
 */
void CollapseExpandFossilBranchProposal::prepareProposal( void )
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
void CollapseExpandFossilBranchProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void CollapseExpandFossilBranchProposal::undoProposal( void )
{
    
    // undo the proposal (only if succeeded)
    if ( failed == false)
    {
        if ( storedNode->isSampledAncestorTip() == true )
        {
            storedNode->setSampledAncestor( false );
            storedNode->getParent().setAge( storedAge );
        }
        else
        {
            storedNode->getParent().setAge( storedAge );
            storedNode->setSampledAncestor( true );
        }
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void CollapseExpandFossilBranchProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tau )
    {
        tau = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if ( oldN == origin )
    {
        origin = static_cast<TypedDagNode<double>* >(newN);
    }
    
}


void CollapseExpandFossilBranchProposal::setProposalTuningParameter(double tp)
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
void CollapseExpandFossilBranchProposal::tune( double rate )
{
    
}

