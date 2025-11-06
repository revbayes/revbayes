/*
 * MarkovTimesDistribution.h
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#ifndef MARKOVTIMESDISTRIBUTION_H_
#define MARKOVTIMESDISTRIBUTION_H_

#include <set>

#include "AbstractEventTimesDistribution.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "OrderedEventTimes.h"

namespace RevBayesCore {

	class MarkovTimesDistribution : public TypedDistribution< OrderedEventTimes >, public AbstractEventTimesDistribution, public MemberObject< std::int64_t >, public MemberObject< RbVector<double> > {

	public:

		MarkovTimesDistribution(const TypedDagNode< double > *rate_, const TypedDagNode< double > *age_);
		MarkovTimesDistribution(const MarkovTimesDistribution &d);
		virtual ~MarkovTimesDistribution();

		MarkovTimesDistribution*                            clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

        // exposed methods
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;             //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const; //!< Map the member methods to internal function calls

        // abstract methods
		void                                                simulateEventTime(double &time, double &ln_prob);
		void                                                getRandomTime(double &time, double &ln_prob);
		void                                                removeTime(double time);
		void                                                addTime(double time);

	protected:

        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter

	private:

        // methods
        OrderedEventTimes*                                  simulate();
        void                                                simulateChildren();

        // members
        const TypedDagNode< double >* rate;
        const TypedDagNode< double >* age;
        size_t                        num_elements;

	};

} /* namespace RevBayesCore */

#endif /* MARKOVTIMESDISTRIBUTION_H_ */
