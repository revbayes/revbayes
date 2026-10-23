#ifndef DagNodeTypeUtilities_H
#define DagNodeTypeUtilities_H

#include "DagNode.h"
#include "TypeUtilities.h"

#include <type_traits>

namespace RevBayesCore {

// Replace a stored DAG-node pointer, preserving null and checking its declared type in debug builds.
// Ownership, graph connections, and any other updates remain the caller's responsibility.
template<class TargetNode, class SourceNode>
    requires std::is_base_of_v<DagNode, std::remove_cv_t<TargetNode>> &&
             std::is_base_of_v<DagNode, std::remove_cv_t<SourceNode>> &&
             requires (SourceNode* replacement) { assume_cast<TargetNode*>(replacement); }
void replaceNodeReference(TargetNode*& node, SourceNode* replacement)
{
    node = assume_cast<TargetNode*>(replacement);
}

}

#endif
