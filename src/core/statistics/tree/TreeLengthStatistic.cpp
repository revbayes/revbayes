#include "TreeLengthStatistic.h"

#include <vector>

#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

TreeLengthStatistic::TreeLengthStatistic(const TypedDagNode<Tree> *t) : TypedFunction<double>( new double(0.0) ), tree( t ) {
    // add the tree parameter as a parent
    addParameter( tree );
    
    update();
}

TreeLengthStatistic::~TreeLengthStatistic( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


TreeLengthStatistic* TreeLengthStatistic::clone( void ) const
{
    return new TreeLengthStatistic( *this );
}

void TreeLengthStatistic::update( void )
{
    
    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
    std::vector<TopologyNode*>::iterator it = nodes.begin();
    std::vector<TopologyNode*>::iterator it_end = nodes.end();
    double treeLength = 0.0;
    for ( ; it != it_end; it++)
    {
        treeLength += tree->getValue().getBranchLengthForNode(*(*it));
    }
    
    *value = treeLength;
}


void TreeLengthStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}
