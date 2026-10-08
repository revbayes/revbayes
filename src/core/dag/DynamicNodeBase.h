#ifndef DynamicNodeBase_H
#define DynamicNodeBase_H

#include <map>
#include <string>

namespace RevBayesCore {

    class DagNode;
    class DagNodeMap;
    class DeterministicNodeBase;
    class StochasticNodeBase;

    /** Implements dynamic-node behavior that does not depend on the node's value type. */
    class DynamicNodeBase {

        friend class DeterministicNodeBase;
        friend class StochasticNodeBase;

    public:
        DynamicNodeBase(void);
        DynamicNodeBase(const DynamicNodeBase &n);
                                                           ~DynamicNodeBase(void) = default;

    protected:
        DagNode*                                           cloneDAG(const DagNode &owner, DagNodeMap &nodesMap, std::map<std::string, const DagNode*> &names) const;
        const std::string&                                 getRevTypeOfValue(void) const;
        bool                                               isTouched(void) const;
        void                                               keepMe(void);
        void                                               restoreMe(void);
        void                                               touchMe(void);

        // members
        bool                                               touched;
    };

}

#endif
