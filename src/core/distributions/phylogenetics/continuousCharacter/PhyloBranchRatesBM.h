#ifndef PhyloBranchRatesBM_H
#define PhyloBranchRatesBM_H

#include "RbVector.h"
#include "TopologyNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class DagNode;
    class Tree;
    
    template <class valueType> class TypedDagNode;
    
    class PhyloBranchRatesBM : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        PhyloBranchRatesBM(const TypedDagNode< Tree > *t, const TypedDagNode< double >* r, const TypedDagNode< double >* s, const TypedDagNode< double >* d);
        
        // public member functions
        PhyloBranchRatesBM*                                     clone(void) const;                                                                      //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        // helper methods
        void                                                    simulate();
        double                                                  recursiveLnProb(const TopologyNode& n, std::vector<double> &parent);
        void                                                    recursiveSimulate(const TopologyNode& n, std::vector<double> &parent);
        
        // private members
        const TypedDagNode< Tree >*                             tau;
        const TypedDagNode< double >*                           root_state;
        const TypedDagNode< double >*                           sigma;
        const TypedDagNode< double >*                           drift;
        
    };
    
}
#endif
