#ifndef SRC_CORE_MONITORS_NEXUSMONITOR_H_
#define SRC_CORE_MONITORS_NEXUSMONITOR_H_

#include <iosfwd>
#include <vector>

#include "AbstractFileMonitor.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;

/** @brief Monitor to output a distribution of trees in Nexus format.
 * The provided branch- or node-specific variables will be output on the trees as metadata.
 *
 * @see ExtendedNewickTreeMonitor for output in table format
* */
class NexusMonitor: public AbstractFileMonitor {
public:
    NexusMonitor(TypedDagNode<Tree> *t, const std::vector<DagNode*> &n, bool np, std::uint64_t g, const std::string &fname,
                 bool ap = false, bool taxa = true);
    NexusMonitor* clone(void) const;

    virtual void printHeader();
    virtual void monitor(std::uint64_t gen);
    virtual void closeStream();    
    void swapNode(DagNode *oldN, DagNode *newN);

protected:
    bool isNodeParameter;  //!< whether data is on the nodes or branches
    bool writeTaxa;  //!< whether to write a taxa block
    TypedDagNode<Tree>* tree;  //!< monitored tree
    std::vector<DagNode*> nodeVariables;  //!< variables associated with the tree
};

} /* namespace RevBayesCore */

#endif /* SRC_CORE_MONITORS_NEXUSMONITOR_H_ */
