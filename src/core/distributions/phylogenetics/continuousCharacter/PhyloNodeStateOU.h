#ifndef PhyloNodeStateOU_H
#define PhyloNodeStateOU_H

#include "RbVector.h"
#include "TopologyNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class DagNode;
    class Tree;
    
    template <class valueType> class TypedDagNode;
    
    class PhyloNodeStateOU : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        PhyloNodeStateOU(const TypedDagNode< Tree > *tr, const TypedDagNode< double >* r, const TypedDagNode< double >* s, const TypedDagNode< double >* a, const TypedDagNode< double >* th);
        
        // public member functions
        PhyloNodeStateOU*                                       clone(void) const;                                                          //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                            getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                            keepSpecialization(const DagNode* affecter);
        virtual void                                            restoreSpecialization(const DagNode *restorer);
        virtual void                                            touchSpecialization(const DagNode *toucher, bool touchAll);

    private:
        // helper methods
        void                                                    simulate();
        double                                                  recursiveLnProb(const TopologyNode& n);
        void                                                    recursiveSimulate(const TopologyNode& n);
        

        // private members
        const TypedDagNode< Tree >*                             tau;
        const TypedDagNode< double >*                           root_state;
        const TypedDagNode< double >*                           sigma;
        const TypedDagNode< double >*                           alpha;
        const TypedDagNode< double >*                           theta;

    };
    
}
#endif
