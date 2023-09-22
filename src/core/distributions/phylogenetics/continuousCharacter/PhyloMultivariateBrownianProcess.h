#ifndef PhyloMultivariateBrownianProcess_H
#define PhyloMultivariateBrownianProcess_H

#include <cstddef>
#include <vector>

#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
class Tree;
    
    class PhyloMultivariateBrownianProcess : public TypedDistribution< RbVector< RbVector<double> > > {
        
    public:
        // constructor(s)
        PhyloMultivariateBrownianProcess(const TypedDagNode< Tree > *intau, const TypedDagNode<MatrixReal>* insigma);
        
        // public member functions
        PhyloMultivariateBrownianProcess*                       clone(void) const;
        
        double                                                  computeLnProbability(void);
        size_t                                                  getDim() const {return sigma->getValue().getDim();}
        void                                                    redrawValue(void);        
        const Tree*                                             getTimeTree() const {return &tau->getValue();}
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

        // special handling of state changes
        void                                                    keepSpecialization(const DagNode* affecter);
        void                                                    restoreSpecialization(const DagNode *restorer);
        void                                                    touchSpecialization(const DagNode *toucher, bool touchAll);
        
    private:
        // helper methods
        void                                                    simulate();
        double                                                  recursiveLnProb(const TopologyNode& n);
        void                                                    recursiveSimulate(const TopologyNode& n);

        // special handling of state changes
        void                                                    flagNodes();        
        void                                                    corruptAll();
        void                                                    recursiveCorruptAll(const TopologyNode& n);
        
        // private members
        const TypedDagNode< Tree >*                             tau;
        const TypedDagNode< MatrixReal >*                       sigma;
        
        std::vector<bool>                                       dirty_nodes;
        std::vector<double>                                     nodeLogProbs;
    };
    
}

#endif /* defined(__revbayes__MultivariateBrownianProcess__) */
