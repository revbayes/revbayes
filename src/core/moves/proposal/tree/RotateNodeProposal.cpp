#include "RotateNodeProposal.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;


RotateNodeProposal::RotateNodeProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    node_index( 0 ),
    failed( false )
{
    addNode( tree );
}


void RotateNodeProposal::cleanProposal( void )
{
}


RotateNodeProposal* RotateNodeProposal::clone( void ) const
{
    return new RotateNodeProposal( *this );
}


const std::string& RotateNodeProposal::getProposalName( void ) const
{
    static std::string name = "RotateNode";

    return name;
}


double RotateNodeProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Rebuild a node's child list in the given order. A node does not expose its child vector, and
 * adding at offset 0 appends, so the children are removed and then re-added in order.
 */
void RotateNodeProposal::setChildren( size_t node_i, const std::vector<TopologyNode*> &c )
{
    TopologyNode &node = tree->getValue().getNode( node_i );

    for (size_t i = 0; i < c.size(); i++)
    {
        node.removeChild( c[i] );
    }
    for (size_t i = 0; i < c.size(); i++)
    {
        node.addChild( c[i], 0 );
    }
}


/**
 * Permute the children of a random internal node. Ages and the clade set are untouched.
 */
double RotateNodeProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    failed = false;

    // which sibling continues the parent's species is the free parameter; child order carries
    // nothing. A sampled ancestor node is skipped, its continuation being forced onto the tip
    std::vector<size_t> candidates;
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        const TopologyNode& n = tau.getNode(i);

        if ( n.isTip() == true || n.getNumberOfChildren() < 2 ) continue;
        if ( n.isSampledAncestorTipOrParent() == true ) continue;

        size_t n_cont = 0;
        for (size_t c = 0; c < n.getNumberOfChildren(); c++)
        {
            if ( n.getChild(c).continuesParentSpecies() == true ) ++n_cont;
        }

        // exactly one holder leaves somewhere to move it to
        if ( n_cont == 1 ) candidates.push_back( i );
    }

    if ( candidates.empty() == true )
    {
        failed = true;

        // reject rather than return 0, or a move that could not act counts as accepted
        return RbConstants::Double::neginf;
    }

    node_index = candidates[ size_t( rng->uniform01() * candidates.size() ) ];

    TopologyNode& node = tau.getNode( node_index );

    std::vector<TopologyNode*> others;
    stored_continuer = NULL;

    for (size_t c = 0; c < node.getNumberOfChildren(); c++)
    {
        TopologyNode *child = &node.getChild(c);

        if ( child->continuesParentSpecies() == true ) stored_continuer = child;
        else                                          others.push_back( child );
    }

    // the reverse move draws from a set of the same size, so the proposal is symmetric
    moved_continuer = others[ size_t( rng->uniform01() * others.size() ) ];

    stored_continuer->setContinuesParentSpecies( false );
    moved_continuer->setContinuesParentSpecies( true );

    return 0.0;
}


void RotateNodeProposal::prepareProposal( void )
{
}


void RotateNodeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void RotateNodeProposal::undoProposal( void )
{
    if ( failed == false )
    {
        moved_continuer->setContinuesParentSpecies( false );
        stored_continuer->setContinuesParentSpecies( true );
    }
}


void RotateNodeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void RotateNodeProposal::setProposalTuningParameter(double tp)
{
}


void RotateNodeProposal::tune( double rate )
{
}
