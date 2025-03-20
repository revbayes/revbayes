#ifndef DebugMove_H
#define DebugMove_H

#include <string>
#include <map>

#include "DagNode.h"
#include "RbOrderedSet.h"

constexpr double rel_err_threshhold = 1.0e-11;
constexpr int err_precision = 11;

void compareNodePrs(const std::string& name, const std::map<const RevBayesCore::DagNode*, double>& pdfs1, const std::map<const RevBayesCore::DagNode*, double>& pdfs2, const std::string& msg);

std::map<const RevBayesCore::DagNode*, double> getNodePrs(const std::vector<RevBayesCore::DagNode*>& nodes, const RevBayesCore::RbOrderedSet<RevBayesCore::DagNode*>& affected_nodes);
#endif
