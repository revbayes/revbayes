/*
 * MarkovTimesDistribution.cpp
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#include "MarkovTimesDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "DistributionExponential.h"
#include "DistributionPoisson.h"
#include "DagNode.h"
#include "StochasticNode.h"

namespace RevBayesCore {

MarkovTimesDistribution::MarkovTimesDistribution(const TypedDagNode< double > *rate_, const TypedDagNode< double > *age_) :
	TypedDistribution< OrderedEventTimes >( new OrderedEventTimes() ),
	rate( rate_ ),
	age( age_ ),
	num_elements(0)
{
	addParameter(rate);
	addParameter(age);
	this->value = simulate();
}

MarkovTimesDistribution::MarkovTimesDistribution(const MarkovTimesDistribution &d) :
	TypedDistribution< OrderedEventTimes >( new OrderedEventTimes() ),
	rate( d.rate ),
	age( d.age ),
	num_elements(d.num_elements)
{
	addParameter(rate);
	addParameter(age);
	this->value = simulate();
}

MarkovTimesDistribution::~MarkovTimesDistribution()
{
}

MarkovTimesDistribution* MarkovTimesDistribution::clone( void ) const
{
    return new MarkovTimesDistribution( *this );
}

double MarkovTimesDistribution::computeLnProbability(void)
{
	// get the number of events
	const size_t &num_events = this->value->size();

	// the number of events is Poisson-distributed
	double ln_prob = RbStatistics::Poisson::lnPdf(rate->getValue() * age->getValue(), num_events);

	return ln_prob;
}

void MarkovTimesDistribution::redrawValue(void)
{
	// simulate the new event times
	this->setValue( simulate() );

	// check whether children need to be simulated
	const std::vector<DagNode*>& children = this->dag_node->getChildren();
	for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		(*it)->redraw();
	}


}





void MarkovTimesDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == rate)
    {
    	rate = static_cast<const TypedDagNode< double > *>( newP );
    }
    else if (oldP == age)
    {
    	age = static_cast<const TypedDagNode< double > *>( newP );
    }
}

OrderedEventTimes* MarkovTimesDistribution::simulate()
{
	// get a random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

	// create the new container
	OrderedEventTimes* ret = new OrderedEventTimes();

	// get the rate and age
	double the_rate = rate->getValue();
	double the_time = age->getValue();

	// simulate forward
	double current_time = 0.0;
	while (true)
	{
		current_time += RbStatistics::Exponential::rv(the_rate, *rng);
		if (current_time > the_time)
		{
			break;
		}
		ret->addEvent(current_time);
	}

	return ret;
}



} /* namespace RevBayesCore */
