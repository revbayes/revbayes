#ifndef PhyloBranchRatesFunction_H
#define PhyloBranchRatesFunction_H

#include <cstddef>
#include <vector>

#include "RbVector.h"
#include "TypedFunction.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class PhyloBranchRatesFunction : public TypedFunction< RbVector<double> > {
        
    public:
        PhyloBranchRatesFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *d, bool as_log);
        virtual                                                ~PhyloBranchRatesFunction(void);                                                     //!< Virtual destructor
        
        // public member functions
        PhyloBranchRatesFunction*                               clone(void) const;                                                                  //!< Create an independent clone
        void                                                    keep(const DagNode* affecter);
        void                                                    restore(const DagNode *restorer);
        void                                                    reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                    touch(const DagNode *toucher );
        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        void                                                    recursiveUpdate(const TopologyNode& n);

        // members
        const TypedDagNode<Tree>*                               tau;
        const TypedDagNode< RbVector<double> >*                 node_state;
        bool                                                    as_log;
        
        std::set<size_t>                                        touched_node_indices;

    };
    
}

#endif

