#ifndef PhyloMultivariateBrownianProcess_H
#define PhyloMultivariateBrownianProcess_H

#include <iostream>


#include "ModelVector.h"
#include "MatrixReal.h"
#include "Real.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
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
        void                                                    keepSpecialization(DagNode* affecter);
        void                                                    restoreSpecialization(DagNode *restorer);
        void                                                    touchSpecialization(DagNode *toucher, bool touchAll);
        
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
