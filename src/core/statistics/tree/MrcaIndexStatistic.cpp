#include "MrcaIndexStatistic.h"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "RbException.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

MrcaIndexStatistic::MrcaIndexStatistic(const TypedDagNode<Tree> *t, const Clade &c) : TypedFunction<std::int64_t>( new std::int64_t(-1) ),
    tree( t ),
    clade( c ),
    index( -1 )
{
    // add the tree parameter as a parent
    addParameter( tree );
    
    initialize();
    update();
}


MrcaIndexStatistic::~MrcaIndexStatistic( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}



MrcaIndexStatistic* MrcaIndexStatistic::clone( void ) const
{
    
    return new MrcaIndexStatistic( *this );
}


void MrcaIndexStatistic::initialize( void )
{
    clade.resetTaxonBitset( tree->getValue().getTaxonBitSetMap() );
    taxa_count = clade.size();
    index = -1;
    
}


void MrcaIndexStatistic::update( void )
{
    if ( not tree->getValue().isRooted() )
    {
        throw RbException("Most recent common ancestor is undefined for unrooted trees.");
    }

    const Tree& t = tree->getValue();
    const TopologyNode& root = t.getRoot();
    const TopologyNode* node = root.getMrca( clade );

    if ( node == NULL )
    {
        throw RbException("MrcaIndex-Statistics can only be applied if clade is present.");
    }

    index = static_cast<int>( node->getIndex() );
    // Return 1-based index so that .getDescendantTaxa() and other Tree methods expecting 1-based node indices work
    *value = index + 1;
}



void MrcaIndexStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
        index = -1;
    }
    
}


