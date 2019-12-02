#ifndef PhyloBrownianProcess_H
#define PhyloBrownianProcess_H

#include "RbVector.h"
#include "Tree.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class PhyloBrownianProcess : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        PhyloBrownianProcess(const TypedDagNode< Tree > *t, const TypedDagNode< double >* s, const TypedDagNode< double >* d);
        
        // public member functions
        PhyloBrownianProcess*                                   clone(void) const;                                                                      //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    redrawValue(void);
        
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
        const TypedDagNode< double >*                           drift;
        
    };
    
}
#endif
