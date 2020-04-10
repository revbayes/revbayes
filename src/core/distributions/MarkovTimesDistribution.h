/*
 * MarkovTimesDistribution.h
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#ifndef MARKOVTIMESDISTRIBUTION_H_
#define MARKOVTIMESDISTRIBUTION_H_

#include <set>

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "OrderedEventTimes.h"

namespace RevBayesCore {

	class MarkovTimesDistribution : public TypedDistribution< OrderedEventTimes > {

	public:

		MarkovTimesDistribution(const TypedDagNode< double > *rate_, const TypedDagNode< double > *age_);
		MarkovTimesDistribution(const MarkovTimesDistribution &d);
		virtual ~MarkovTimesDistribution();

		MarkovTimesDistribution*                            clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

	protected:

        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter

	private:

        // methods
        OrderedEventTimes*                                  simulate();

        // members
        const TypedDagNode< double >* rate;
        const TypedDagNode< double >* age;
        size_t                        num_elements;

	};

} /* namespace RevBayesCore */

#endif /* MARKOVTIMESDISTRIBUTION_H_ */
