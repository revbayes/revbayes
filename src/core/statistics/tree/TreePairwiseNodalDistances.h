//
//  TreePairwiseNodalDistances.h
//  Will Freyman 5/21/2018
//

#ifndef TreePairwiseNodalDistances_h
#define TreePairwiseNodalDistances_h

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "DistanceMatrix.h"

namespace RevBayesCore {
class DagNode;
class Tree;
    
	class TreePairwiseNodalDistances : public TypedFunction<RevBayesCore::DistanceMatrix> {
        
    public:
        TreePairwiseNodalDistances(const TypedDagNode<Tree> *t);                                                                            
        virtual                                    ~TreePairwiseNodalDistances(void);                                                       
        
        TreePairwiseNodalDistances&                 operator=(const TreePairwiseNodalDistances& t);
        TreePairwiseNodalDistances*                 clone(void) const;                                                                      
        void                                        update(void);                                                                           
        
    protected:

        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        
        
    private:
		
        RevBayesCore::DistanceMatrix* getDistanceMatrix(const TypedDagNode<Tree>& tree);
        const TypedDagNode<Tree>*                   tree;
        
    };
    
}


#endif 
