#ifndef DeterministicNodeBase_H
#define DeterministicNodeBase_H

#include "RbOrderedSet.h"

#include <vector>

namespace RevBayesCore {

    class DagNode;
    class DynamicNodeBase;
    class Function;

    /** Implements deterministic-node behavior that does not depend on the node's value type. */
    class DeterministicNodeBase {

    public:
        explicit DeterministicNodeBase(Function *f);
        DeterministicNodeBase(const DeterministicNodeBase &n);
                                                           ~DeterministicNodeBase(void);

    protected:
        void                                               assign(const DeterministicNodeBase &n, DagNode &owner);
        void                                               attachToFunctionParameters(DagNode &owner) const;
        void                                               detachFromFunctionParameters(DagNode &owner) const;
        void                                               getAffected(DagNode &owner, RbOrderedSet<DagNode*> &affected) const;
        void                                               getIntegratedParents(const DagNode &owner, RbOrderedSet<DagNode*> &integratedParents) const;
        std::vector<const DagNode*>                        getParents(void) const;
        bool                                               isConstant(void) const;
        void                                               keepMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *affecter);
        void                                               reInitializeMe(void);
        void                                               restoreMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *restorer);
        void                                               swapParent(DagNode &owner, const DagNode *oldParent, const DagNode *newParent);
        void                                               touchMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *toucher, bool touchAll);

        // members
        Function                                          *function;
        mutable bool                                       needs_update;
        bool                                               force_update;
    };

}

#endif
