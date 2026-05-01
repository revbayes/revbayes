#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <optional>
#include <vector>
#include "DagNode.h"

namespace Transforms
{
    std::optional<double> add_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> add_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_add_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> mul_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> mul_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_mul_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> exp_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> exp_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_exp_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> exp_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> exp_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_exp_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> log_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_log_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> logit_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> logit_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_logit_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> sub1_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> sub1_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_sub1_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> sub2_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> sub2_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_sub2_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

    std::optional<double> invlogit_transform(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> invlogit_inverse(const std::vector<const RevBayesCore::DagNode*>& params, double x);
    std::optional<double> log_invlogit_prime(const std::vector<const RevBayesCore::DagNode*>& params, double x);

};

#endif
