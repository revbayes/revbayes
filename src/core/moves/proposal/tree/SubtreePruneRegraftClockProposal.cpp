#include "SubtreePruneRegraftClockProposal.h"

#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;


SubtreePruneRegraftClockProposal::SubtreePruneRegraftClockProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    stored_node( NULL ),
    stored_sibling( NULL ),
    stored_age( 0.0 ),
    stored_slot( 0 ),
    failed( false )
{
    addNode( tree );
}


void SubtreePruneRegraftClockProposal::cleanProposal( void )
{
}


SubtreePruneRegraftClockProposal* SubtreePruneRegraftClockProposal::clone( void ) const
{
    return new SubtreePruneRegraftClockProposal( *this );
}


const std::string& SubtreePruneRegraftClockProposal::getProposalName( void ) const
{
    static std::string name = "SubtreePruneRegraftClock";

    return name;
}


double SubtreePruneRegraftClockProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


void SubtreePruneRegraftClockProposal::markSubtree( const TopologyNode &n, std::vector<bool> &in ) const
{
    in[ n.getIndex() ] = true;

    for (size_t i = 0; i < n.getNumberOfChildren(); i++)
    {
        markSubtree( n.getChild(i), in );
    }
}


/** A prunable node has a parent to detach from and a grandparent to reconnect its sibling to. */
size_t SubtreePruneRegraftClockProposal::countMovable( void ) const
{
    const Tree &tau = tree->getValue();

    size_t n = 0;
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        const TopologyNode &nd = tau.getNode(i);

        if ( nd.isRoot() == false && nd.getParent().isRoot() == false ) n++;
    }

    return n;
}


/** The ages the branch above target can give to a node joining it to sub. */
double SubtreePruneRegraftClockProposal::window( const TopologyNode &sub, const TopologyNode &target ) const
{
    double lo = target.getAge() > sub.getAge() ? target.getAge() : sub.getAge();
    double hi = target.getParent().getAge();

    return hi - lo;
}


/** Branches that can take the subtree: not inside it, not its own stem, and old enough. */
void SubtreePruneRegraftClockProposal::attachments( const TopologyNode &sub, std::vector<TopologyNode*> &out ) const
{
    Tree &tau = tree->getValue();

    std::vector<bool> inside( tau.getNumberOfNodes(), false );
    markSubtree( sub, inside );

    const TopologyNode &parent = sub.getParent();

    out.clear();
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        TopologyNode &nd = tau.getNode(i);

        if ( nd.isRoot() == true || inside[i] == true ) continue;
        if ( &nd == &parent ) continue;

        if ( window( sub, nd ) > 0.0 ) out.push_back( &nd );
    }
}


/** Detach sub with its parent, closing the gap behind it. Each node keeps its child slot. */
void SubtreePruneRegraftClockProposal::prune( TopologyNode *sub )
{
    TopologyNode *parent  = &sub->getParent();
    TopologyNode *gparent = &parent->getParent();

    TopologyNode *sibling = &parent->getChild( 0 );
    if ( sibling == sub ) sibling = &parent->getChild( 1 );

    size_t pos = gparent->removeChild( parent );
    stored_slot = parent->removeChild( sibling );
    gparent->addChild( sibling, pos );
    sibling->setParent( gparent );
}


/** Split the branch above target at age and hang the pruned subtree there. */
void SubtreePruneRegraftClockProposal::attach( TopologyNode *sub, TopologyNode *target, double age )
{
    TopologyNode *parent  = &sub->getParent();
    TopologyNode *gparent = &target->getParent();

    size_t pos = gparent->removeChild( target );
    gparent->addChild( parent, pos );
    parent->setParent( gparent );
    parent->addChild( target, stored_slot );
    target->setParent( parent );

    parent->setAge( age );
}


double SubtreePruneRegraftClockProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    failed = false;

    size_t movable_before = countMovable();
    if ( movable_before == 0 )
    {
        failed = true;

        return 0.0;
    }

    // pick the subtree to move
    TopologyNode *node = NULL;
    size_t seen = 0;
    size_t pick = size_t( rng->uniform01() * movable_before );
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        TopologyNode &nd = tau.getNode(i);

        if ( nd.isRoot() == true || nd.getParent().isRoot() == true ) continue;
        if ( seen == pick ) { node = &nd; break; }
        seen++;
    }
    if ( node == NULL )
    {
        failed = true;

        return 0.0;
    }

    stored_node = node;
    stored_sibling = &node->getParent().getChild( 0 );
    if ( stored_sibling == node ) stored_sibling = &node->getParent().getChild( 1 );
    stored_age = node->getParent().getAge();

    // pruning the subtree out of either state gives the same tree, so the attachment set and
    // both age windows are measured on it and the combinatorial term cancels
    prune( node );

    std::vector<TopologyNode*> targets;
    attachments( *node, targets );

    double rev_window = window( *node, *stored_sibling );

    if ( targets.empty() == true || rev_window <= 0.0 )
    {
        attach( node, stored_sibling, stored_age );
        failed = true;

        return 0.0;
    }

    TopologyNode *target = targets[ size_t( rng->uniform01() * targets.size() ) ];

    double lo = target->getAge() > node->getAge() ? target->getAge() : node->getAge();
    double fwd_window = window( *node, *target );
    double new_age = lo + rng->uniform01() * fwd_window;

    attach( node, target, new_age );

    size_t movable_after = countMovable();
    if ( movable_after == 0 )
    {
        return RbConstants::Double::neginf;
    }

    return log( double(movable_before) / double(movable_after) )
         + log( fwd_window / rev_window );
}


void SubtreePruneRegraftClockProposal::prepareProposal( void )
{
}


void SubtreePruneRegraftClockProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void SubtreePruneRegraftClockProposal::undoProposal( void )
{
    if ( failed == false )
    {
        prune( stored_node );
        attach( stored_node, stored_sibling, stored_age );
    }
}


void SubtreePruneRegraftClockProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void SubtreePruneRegraftClockProposal::setProposalTuningParameter(double tp)
{
}


void SubtreePruneRegraftClockProposal::tune( double rate )
{
}
