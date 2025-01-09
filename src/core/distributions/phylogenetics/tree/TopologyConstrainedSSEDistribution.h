/*
 * TopologyConstrainedSSEDistribution.h
 *
 *  Created on: Jan 31, 2024
 *      Author: mike
 */

#ifndef TOPOLOGYCONSTRAINEDSSEDISTRIBUTION_H_
#define TOPOLOGYCONSTRAINEDSSEDISTRIBUTION_H_

#include "RbVector.h"
#include "TopologyConstrainedTreeDistribution.h"

namespace RevBayesCore {

class TopologyConstrainedSSEDistribution: public TopologyConstrainedTreeDistribution, public MemberObject< RbVector<double> > {

public:

    TopologyConstrainedSSEDistribution(TypedDistribution<Tree>* base_dist, const std::vector<Clade> &c, Tree *t, long age_check_precision);
	TopologyConstrainedSSEDistribution(const TopologyConstrainedSSEDistribution &d);

	virtual ~TopologyConstrainedSSEDistribution();

    void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;
    RevLanguage::RevPtr<RevLanguage::RevVariable>       executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);


};

} /* namespace RevBayesCore */

#endif /* TOPOLOGYCONSTRAINEDSSEDISTRIBUTION_H_ */
