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

    std::vector<size_t> internal;
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        const TopologyNode& n = tau.getNode(i);

        if ( n.isTip() == false && n.getNumberOfChildren() > 1 )
        {
            internal.push_back( i );
        }
    }

    if ( internal.empty() == true )
    {
        failed = true;

        return 0.0;
    }

    node_index = internal[ size_t( rng->uniform01() * internal.size() ) ];

    stored_children = tau.getNode( node_index ).getChildren();

    // an order drawn uniformly from those that differ from the current one; the reverse draw has
    // the same probability, so the proposal is symmetric
    std::vector<TopologyNode*> rotated = stored_children;
    for (size_t attempt = 0; attempt < 100; attempt++)
    {
        for (size_t i = rotated.size() - 1; i > 0; i--)
        {
            size_t j = size_t( rng->uniform01() * (i + 1) );
            if ( j > i ) j = i;

            TopologyNode *tmp = rotated[i];
            rotated[i] = rotated[j];
            rotated[j] = tmp;
        }
        if ( rotated != stored_children ) break;
    }

    if ( rotated == stored_children )
    {
        failed = true;

        return 0.0;
    }

    setChildren( node_index, rotated );

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
        setChildren( node_index, stored_children );
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
