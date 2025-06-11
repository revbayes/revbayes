#ifndef DebugMove_H
#define DebugMove_H

#include <string>
#include <map>

#include "DagNode.h"
#include "RbOrderedSet.h"
#include <variant>
#include <iostream>

constexpr double rel_err_threshhold = 1.0e-11;
constexpr int err_precision = 11;

struct ProbOrError: public std::variant<double, std::string>
{
    const double* is_prob() const;
    const std::string* is_error() const;
    double is_prob_or(double d) const;

    using variant::variant;
    ProbOrError() = delete;
};

std::ostream& operator<<(std::ostream& o, const ProbOrError& e);

typedef std::map<const RevBayesCore::DagNode*, ProbOrError> NodePrMap;

void compareNodePrs(const std::string& name, const NodePrMap& pdfs1, const NodePrMap& pdfs2, const std::string& msg);

NodePrMap getNodePrs(const std::vector<RevBayesCore::DagNode*>& nodes, const RevBayesCore::RbOrderedSet<RevBayesCore::DagNode*>& affected_nodes);
#endif
