#include "TreePairwiseDistances.h"

#include "Tree.h"
#include "TreeUtilities.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

TreePairwiseDistances::TreePairwiseDistances(const TypedDagNode<Tree> *t) : TypedFunction<RevBayesCore::DistanceMatrix>( new RevBayesCore::DistanceMatrix(t->getValue().getNumberOfTips()) ),
    tree( t )
{
    // add the tree parameter as a parent
    addParameter( tree );
    
    update();
}

TreePairwiseDistances::~TreePairwiseDistances( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


TreePairwiseDistances* TreePairwiseDistances::clone( void ) const
{
    return new TreePairwiseDistances( *this );
}



RevBayesCore::DistanceMatrix* TreePairwiseDistances::getDistanceMatrix(const TypedDagNode<Tree>& tree)
{
    
    RevBayesCore::DistanceMatrix* matrix = TreeUtilities::getDistanceMatrix ( tree.getValue() );
    
    /*new RevBayesCore::MatrixReal( tree.getValue().getNumberOfTips() );
     
     std::vector<std::string> names = tree.getValue().getTipNames( ) ;
     
     
     
     std::map< std::string, int > namesToId;
     
     for (size_t i = 0; i < names.size(); ++i) {
     namesToId[ names[i] ] = i;
     }
     
     std::vector< std::pair<std::string, double> > distsToRoot;
     
     processDistsInSubtree( tree.getValue().getRoot() , *matrix, distsToRoot, namesToId);*/
    
    return matrix;
}




void TreePairwiseDistances::update( void )
{
    
    delete value;
    value = getDistanceMatrix( *tree );
}


void TreePairwiseDistances::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}
