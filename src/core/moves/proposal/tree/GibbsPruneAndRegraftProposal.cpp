#include "GibbsPruneAndRegraftProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
GibbsPruneAndRegraftProposal::GibbsPruneAndRegraftProposal( StochasticNode<Tree> *n ) : Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void GibbsPruneAndRegraftProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
GibbsPruneAndRegraftProposal* GibbsPruneAndRegraftProposal::clone( void ) const
{
    
    return new GibbsPruneAndRegraftProposal( *this );
}



void GibbsPruneAndRegraftProposal::findNewBrothers(std::vector<TopologyNode *> &b, TopologyNode &p, TopologyNode *n)
{
    // security check that I'm not a tip
    if (!n->isTip() && &p != n)
    {
        // check the first child
        TopologyNode* child = &n->getChild( 0 );
        if ( child->getAge() < p.getAge() )
        {
            b.push_back( child );
        }
        else
        {
            findNewBrothers(b, p, child);
        }
        
        // check the second child
        child = &n->getChild( 1 );
        if ( child->getAge() < p.getAge() )
        {
            b.push_back( child );
        }
        else
        {
            findNewBrothers(b, p, child);
        }
        
    }
    
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& GibbsPruneAndRegraftProposal::getProposalName( void ) const
{
    static std::string name = "GPR";
    
    return name;
}


double GibbsPruneAndRegraftProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double GibbsPruneAndRegraftProposal::doProposal( void )
{
    
    // reset flags
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    if ( tau.getNumberOfTips() < 3)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // potential affected nodes for likelihood computation
    RbOrderedSet<DagNode *> affected;
    variable->initiateGetAffectedNodes( affected );
    
    double backward_likelihood = variable->getLnProbability();
    for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
    {
        backward_likelihood += (*it)->getLnProbability();
    }
    int offset = (int) -backward_likelihood;
    double backward = exp(backward_likelihood + offset);
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    TopologyNode* parent        = &node->getParent();
    TopologyNode& grandparent   = parent->getParent();
    TopologyNode& brother       = parent->getChild( 0 );
    // check if we got the correct child
    if ( &brother == node )
    {
        brother = parent->getChild( 1 );
    }
    
    // collect the possible reattachement points
    std::vector<TopologyNode*> new_brothers;
    findNewBrothers(new_brothers, *parent, &tau.getRoot());
    
    // we only need to propose a new tree if there are any other re-attachement points
    if ( new_brothers.size() < 1)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    std::vector<double> weights = std::vector<double>(new_brothers.size(), 0.0);
    double sumOfWeights = 0.0;
    for (size_t i = 0; i<new_brothers.size(); ++i)
    {
        // get the new brother
        TopologyNode* newBro = new_brothers[i];
        
        // do the proposal
        TopologyNode *newGrandparent = pruneAndRegraft(&brother, newBro, parent, grandparent);
        
        // flag for likelihood recomputation
        variable->touch();
        
        // compute the likelihood of the new value
        double priorRatio = variable->getLnProbability();
        double likelihoodRatio = 0.0;
        for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
        {
            likelihoodRatio += (*it)->getLnProbability();
        }
        weights[i] = exp(priorRatio + likelihoodRatio + offset);
        sumOfWeights += weights[i];
        
        // undo proposal
        pruneAndRegraft(newBro, &brother, parent, *newGrandparent);
        
        // restore the previous likelihoods;
        variable->restore();
    }
    
    if (sumOfWeights <= 1E-100 || sumOfWeights != sumOfWeights)
    {
        // hack
        // the proposals have such a small likelihood that they can be neglected
        // or sumOfWeights is NaN 
        //        throw new OperatorFailedException("Couldn't find another proposal with a decent likelihood.");
        failed = true;
        
        return 0.0;
    }
    
    double ran = rng->uniform01() * sumOfWeights;
    size_t index = 0;
    while (ran > 0.0)
    {
        ran -= weights[index];
        index++;
    }
    index--;
    
    TopologyNode* newBro = new_brothers[index];
    
    // now we store all necessary values
    storedBrother       = &brother;
    storedNewBrother    = newBro;
    
    pruneAndRegraft(&brother, newBro, parent, grandparent);
    
    double forward = weights[index];
    
    double forwardProb = (forward / sumOfWeights);
    double backwardProb = (backward / (sumOfWeights - forward + backward));
    double hastingsRatio = log(backwardProb / forwardProb);
    
    return hastingsRatio;
}


/**
 *
 */
void GibbsPruneAndRegraftProposal::prepareProposal( void )
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
void GibbsPruneAndRegraftProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


TopologyNode* GibbsPruneAndRegraftProposal::pruneAndRegraft(TopologyNode *brother, TopologyNode *newBrother, TopologyNode *parent, TopologyNode &grandparent)
{
    
    // prune
    grandparent.removeChild( parent );
    parent->removeChild( brother );
    grandparent.addChild( brother );
    brother->setParent( &grandparent );
    
    // regraft
    TopologyNode* newGrandParent = &newBrother->getParent();
    newGrandParent->removeChild( newBrother );
    newGrandParent->addChild( parent );
    parent->setParent( newGrandParent );
    parent->addChild( newBrother );
    newBrother->setParent( parent );
    
    return newGrandParent;
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void GibbsPruneAndRegraftProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( !failed )
    {
        // undo the proposal
        TopologyNode& parent = storedNewBrother->getParent();
        TopologyNode& newGrandparent = parent.getParent();
        TopologyNode& grandparent = storedBrother->getParent();
        
        // prune
        newGrandparent.removeChild( &parent );
        parent.removeChild( storedNewBrother );
        newGrandparent.addChild( storedNewBrother );
        storedNewBrother->setParent( &newGrandparent );
        
        
        // regraft
        grandparent.removeChild( storedBrother );
        parent.addChild( storedBrother );
        storedBrother->setParent( &parent );
        grandparent.addChild( &parent );
        parent.setParent( &grandparent );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GibbsPruneAndRegraftProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void GibbsPruneAndRegraftProposal::setProposalTuningParameter(double tp)
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
void GibbsPruneAndRegraftProposal::tune( double rate )
{
    
    // nothing to tune
    
}

