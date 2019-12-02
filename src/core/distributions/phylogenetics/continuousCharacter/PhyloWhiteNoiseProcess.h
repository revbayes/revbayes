#ifndef PhyloWhiteNoiseProcess_H
#define PhyloWhiteNoiseProcess_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class PhyloWhiteNoiseProcess : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        PhyloWhiteNoiseProcess(const TypedDagNode< Tree > *t, const TypedDagNode< double >* s);
        
        // public member functions
        PhyloWhiteNoiseProcess*                                 clone(void) const;                                                                      //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    redrawValue(void);
 
        /*
        // special handling of state changes
        void                                                    getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);                          //!< get affected nodes
        void                                                    keepSpecialization(DagNode* affecter);
        void                                                    restoreSpecialization(DagNode *restorer);
        void                                                    touchSpecialization(DagNode *toucher, bool touchAll);
        */
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        // helper methods
        void                                                    simulate();
        double                                                  recursiveLnProb(const TopologyNode& n);
        void                                                    recursiveSimulate(const TopologyNode& n);
        
        // private members
        const TypedDagNode< Tree >*                             tau;
        const TypedDagNode< double >*                           sigma;
        
    };
    
}
#endif
