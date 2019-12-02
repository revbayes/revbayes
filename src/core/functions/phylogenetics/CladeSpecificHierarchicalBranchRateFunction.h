#ifndef CladeSpecificHierarchicalBranchRateFunction_H
#define CladeSpecificHierarchicalBranchRateFunction_H

#include "Clade.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class CladeSpecificHierarchicalBranchRateFunction : public TypedFunction< RbVector<double> > {
        
    public:
        CladeSpecificHierarchicalBranchRateFunction(const TypedDagNode<Tree> *t, const std::vector<Clade> &c, const TypedDagNode< RbVector<double> > *br, const TypedDagNode< RbVector<double> > *cr);
        virtual                                                ~CladeSpecificHierarchicalBranchRateFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        CladeSpecificHierarchicalBranchRateFunction*            clone(void) const;                                                                  //!< Create an independent clone
        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                    flattenClades(void);
        
    private:
        // members
        const TypedDagNode<Tree>*                               tau;
        std::vector<Clade>                                      clades;
        const TypedDagNode< RbVector<double> > *                base_rates;
        const TypedDagNode< RbVector<double> > *                clade_rates;
        size_t                                                  n_branches;
        size_t                                                  n_clades;
        std::vector<size_t>                                     flat_clades;
    };
    
}

#endif

