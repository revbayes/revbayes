#include "SubtreePruneRegraftProposal.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SubtreePruneRegraftProposal::SubtreePruneRegraftProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreePruneRegraftProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreePruneRegraftProposal* SubtreePruneRegraftProposal::clone( void ) const
{
    
    return new SubtreePruneRegraftProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreePruneRegraftProposal::getProposalName( void ) const
{
    static std::string name = "SubtreePruneRegraft";
    
    return name;
}


double SubtreePruneRegraftProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


bool SubtreePruneRegraftProposal::isDescendant(const TopologyNode &n, const TopologyNode &p) {

    if ( n.isRoot() ) {
        return false;
    }

    if ( &n == &p ) {
        return true;
    }

    return isDescendant(n.getParent(), p);
}


/**
 * Perform the proposal.
 *
 * A SPR proposal.
 *
 * \return The hastings ratio.
 */
double SubtreePruneRegraftProposal::doProposal( void )
{
    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // when the tree has fewer than 4 tips, nothing should be done
    size_t num_tips = tau.getNumberOfTips();
    size_t num_root_children = tau.getRoot().getNumberOfChildren();
    if ( num_tips <= num_root_children)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    // now we store all necessary values
    storedChoosenNode   = node;
    TopologyNode &parent = node->getParent();
    TopologyNode &grandparent = parent.getParent();
    storedBrother = &parent.getChild( 0 );
    
    // check if we got the correct child
    if ( node == storedBrother ) {
        storedBrother = &parent.getChild( 1 );
    }
    
    // pick a random new parent node
    TopologyNode* newBrother;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        newBrother = &tau.getNode(index);
    } while ( newBrother->isRoot() || newBrother == &parent || newBrother == storedBrother || newBrother == storedChoosenNode || &(newBrother->getParent()) == storedChoosenNode);
    
    TopologyNode &newGrandparent = newBrother->getParent();
    
    // if the regrafting subtree contains the root
    if (!isDescendant(*newBrother,*node))
    {
        prunedroot = false;
        // now prune
        grandparent.removeChild( &parent );
        parent.removeChild( storedBrother );
        grandparent.addChild( storedBrother );
        storedBrother->setParent( &grandparent, false );

        // re-attach
        newGrandparent.removeChild( newBrother );
        parent.addChild( newBrother );
        newGrandparent.addChild( &parent );
        parent.setParent( &newGrandparent, false );
        newBrother->setParent( &parent, false );
    }
    // if the pruned subtree contains the root
    else
    {
        prunedroot = true;
        // prune
        std::vector<TopologyNode*> storedchildren = storedChoosenNode->getChildren();
        std::vector<size_t> storedindices;
        std::vector<double> storedbrlens;
        for (size_t i = 0 ; i < storedchildren.size(); i++)
        {
            storedindices.push_back(storedchildren[i]->getIndex());
            storedbrlens.push_back(storedchildren[i]->getBranchLength());
            storedChoosenNode->removeChild(storedchildren[i]);
            storedchildren[i]->setParent(NULL, false);
        }
        TopologyNode* newParent = &(newBrother->getParent());
        TopologyNode& regraft = tau.reverseParentChild(*newParent);

        for (size_t i = 0 ; i < storedchildren.size(); i++)
        {
            if (&regraft != storedchildren[i])
            {
                regraft.addChild(storedchildren[i]);
                storedchildren[i]->setParent(&regraft, false);

                storedBrother = storedchildren[i];
            }
            else
            {
                newParent->setIndex(storedindices[i]);
                newParent->setBranchLength(storedbrlens[i]);
            }
        }

        newParent->removeChild(newBrother);

        // re-attach
        storedChoosenNode->addChild(newBrother);
        newBrother->setParent(storedChoosenNode, false);
        storedChoosenNode->addChild(newParent);
        newParent->setParent(storedChoosenNode, false);

        tau.orderNodesByIndex();
    }
    
    return 0.0;
    
}


/**
 *
 */
void SubtreePruneRegraftProposal::prepareProposal( void )
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
void SubtreePruneRegraftProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SubtreePruneRegraftProposal::undoProposal( void )
{

    // we undo the proposal only if it didn't fail
    if ( failed == true )
    {
        return;
    }
    
    // undo the proposal
    if (!prunedroot)
    {
        TopologyNode &parent = storedChoosenNode->getParent();
        TopologyNode &grandparent = parent.getParent();
        TopologyNode* oldBrother = &parent.getChild( 0 );
        TopologyNode &newGrandparent = storedBrother->getParent();

        // check if we got the correct child
        if ( storedChoosenNode == oldBrother ) {
            oldBrother = &parent.getChild( 1 );
        }

        // now prune
        grandparent.removeChild( &parent );
        parent.removeChild( oldBrother );
        grandparent.addChild( oldBrother );
        oldBrother->setParent( &grandparent, false );

        // re-attach
        newGrandparent.removeChild( storedBrother );
        parent.addChild( storedBrother );
        newGrandparent.addChild( &parent );
        parent.setParent( &newGrandparent, false );
        storedBrother->setParent( &parent, false );
    }
    else
    {
        Tree& tau = tree->getValue();

        // prune
        std::vector<TopologyNode*> storedchildren = storedChoosenNode->getChildren();
        std::vector<size_t> storedindices;
        std::vector<double> storedbrlens;
        for (size_t i = 0 ; i < storedchildren.size(); i++)
        {
            storedindices.push_back(storedchildren[i]->getIndex());
            storedbrlens.push_back(storedchildren[i]->getBranchLength());
            storedChoosenNode->removeChild(storedchildren[i]);
            storedchildren[i]->setParent(NULL, false);
        }

        TopologyNode* newParent = &(storedBrother->getParent());
        TopologyNode& oldBrother = tau.reverseParentChild(*newParent);

        for (size_t i = 0 ; i < storedchildren.size(); i++)
        {
            if (&oldBrother != storedchildren[i])
            {
                oldBrother.addChild(storedchildren[i]);
                storedchildren[i]->setParent(&oldBrother, false);
            }
            else
            {
                newParent->setIndex(storedindices[i]);
                newParent->setBranchLength(storedbrlens[i]);
            }
        }

        newParent->removeChild(storedBrother);

        // re-attach
        storedChoosenNode->addChild(storedBrother);
        storedBrother->setParent(storedChoosenNode, false);
        storedChoosenNode->addChild(newParent);
        newParent->setParent(storedChoosenNode, false);

        tau.orderNodesByIndex();
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreePruneRegraftProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void SubtreePruneRegraftProposal::setProposalTuningParameter(double tp)
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
void SubtreePruneRegraftProposal::tune( double rate )
{
    
    // nothing to tune
    
}

