/*
 * TopologyConstrainedSSEDistribution.cpp
 *
 *  Created on: Jan 31, 2024
 *      Author: mike
 */

#include "TopologyConstrainedSSEDistribution.h"

namespace RevBayesCore {

TopologyConstrainedSSEDistribution::TopologyConstrainedSSEDistribution(TypedDistribution<Tree>* base_dist, const std::vector<Clade> &c, Tree *t, long age_check_precision) :
    TopologyConstrainedTreeDistribution(base_dist, c, t, age_check_precision)
{
}

// copy constructor needs to call base copy constructor
TopologyConstrainedSSEDistribution::TopologyConstrainedSSEDistribution(const TopologyConstrainedSSEDistribution &d) : TopologyConstrainedTreeDistribution(d)
{
}

TopologyConstrainedSSEDistribution::~TopologyConstrainedSSEDistribution()
{
}

void TopologyConstrainedSSEDistribution::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
{
}

RevLanguage::RevPtr<RevLanguage::RevVariable> TopologyConstrainedSSEDistribution::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
{
	return base_distribution->executeProcedure(name, args, found);
}

} /* namespace RevBayesCore */
