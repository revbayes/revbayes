#include "TreePairwiseNodalDistances.h"


using namespace RevBayesCore;


TreePairwiseNodalDistances::TreePairwiseNodalDistances(const TypedDagNode<Tree> *t) : TypedFunction<RevBayesCore::DistanceMatrix>( new RevBayesCore::DistanceMatrix(t->getValue().getNumberOfTips()) ),
    tree( t )
{
    addParameter( tree );
    
    update();
}

TreePairwiseNodalDistances::~TreePairwiseNodalDistances( void )
{
}


TreePairwiseNodalDistances* TreePairwiseNodalDistances::clone( void ) const
{
    return new TreePairwiseNodalDistances( *this );
}


RevBayesCore::DistanceMatrix* TreePairwiseNodalDistances::getDistanceMatrix(const TypedDagNode<Tree>& tree)
{
    
    RevBayesCore::DistanceMatrix* matrix = TreeUtilities::getNodalDistanceMatrix ( tree.getValue() );
    
    return matrix;
}




void TreePairwiseNodalDistances::update( void )
{
    
    RevBayesCore::DistanceMatrix matrix =	*(getDistanceMatrix( *tree ) );
    *value = matrix;
}


void TreePairwiseNodalDistances::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}
