#ifndef DebugMove_H
#define DebugMove_H

#include <string>
#include <map>

#include "DagNode.h"
#include "RbOrderedSet.h"
#include <variant>
#include <iostream>

#ifndef NDEBUG
constexpr bool debug_build = true;
#else
constexpr bool debug_build = false;
#endif

constexpr double rel_err_threshhold = 1.0e-11;
constexpr int err_precision = 11;

struct ProbOrError: public std::variant<LogDensity, std::string>
{
    const LogDensity* is_prob() const;
    const std::string* is_error() const;
    LogDensity is_prob_or(LogDensity d) const;
    LogDensity is_prob_or_nan() const;

    using variant::variant;
    ProbOrError() = delete;
};

std::ostream& operator<<(std::ostream& o, const ProbOrError& e);

typedef std::map<const RevBayesCore::DagNode*, ProbOrError> NodePrMap;

void compareNodePrs(const std::string& name, const NodePrMap& pdfs1, const NodePrMap& pdfs2, const std::string& msg);

NodePrMap getNodePrs(const std::vector<RevBayesCore::DagNode*>& nodes, const RevBayesCore::RbOrderedSet<RevBayesCore::DagNode*>& affected_nodes);

void showNodePrs(const std::string& stage, const NodePrMap& pdfs);
void showNodePrs(const std::string& stage, const NodePrMap& pdfs, const std::vector<const RevBayesCore::DagNode*>& nodes);

void checkNotTouched(const std::vector<RevBayesCore::DagNode*>& nodes, const RevBayesCore::RbOrderedSet<RevBayesCore::DagNode*>& affected_nodes);
#endif
