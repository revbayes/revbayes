/*
 * MarkovTimesDistribution.cpp
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#include "MarkovTimesDistribution.h"

#include <set>

#include "AbstractEventsDistribution.h"
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

	// make sure the events are between 0 and age
	const double& max = age->getValue();
	double min = 0.0;
	if ( num_events > 0 )
	{
		const std::set<double>& event_times = this->value->getEventTimes();
		std::set<double>::const_iterator first = event_times.begin();
		std::set<double>::const_iterator last  = first;
		std::advance(last, num_events - 1);

		if( *first < min || *last > max )
		{
			return RbConstants::Double::neginf;
		}
	}

	// the number of events is Poisson-distributed
	double ln_prob = RbStatistics::Poisson::lnPdf(rate->getValue() * age->getValue(), num_events);

	// the time of events is uniformly distributed
	ln_prob -= double(num_events) * std::log(max);

	return ln_prob;
}

void MarkovTimesDistribution::redrawValue(void)
{
	// simulate the new event times
	this->setValue( simulate() );

	// simulate children
	this->simulateChildren();
}

void MarkovTimesDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::int64_t &rv) const
{
    if ( n == "getNumberOfEvents" )
    {
        rv = value->size();
    }
    else
    {
        throw RbException("The markov event model does not have a member method called '" + n + "'.");
    }
}

void MarkovTimesDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    if ( n == "getTimes" )
    {
    	// get the times as a set
    	const std::set<double>& times = value->getEventTimes();

    	// convert the times to a vector
    	std::vector<double> res;
    	for( std::set<double>::const_iterator it = times.begin(); it != times.end(); ++it)
    	{
    		res.push_back( *it );
    	}

    	// add the times to the rv
    	rv = res;
    }
    else
    {
        throw RbException("The markov event model does not have a member method called '" + n + "'.");
    }
}

void MarkovTimesDistribution::simulateEventTime(double &time, double &ln_prob)
{
	// get a random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // compute the boundaries
	const double& max = age->getValue();

	// generate the new age
	double new_age = rng->uniform01() * max;

	// store the values
	time = new_age;

	// store the probability
	ln_prob = -std::log(max);
}

void MarkovTimesDistribution::getRandomTime(double &time, double &ln_prob)
{
	// get a random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the event times
	const std::set<double>& event_times = this->value->getEventTimes();

    // choose an element to update
    size_t num_events = this->value->size();
	size_t index = size_t(rng->uniform01() * num_events);

	// get the value of the chosen element
	std::set<double>::const_iterator the_event = event_times.begin();
	std::advance(the_event, index);

	// store the time
	time = *the_event;

    // compute the boundaries
	const double& max = age->getValue();

	// store the probability
	ln_prob = -std::log(max);
}

void MarkovTimesDistribution::removeTime(double time)
{
	// remove the event
	bool success = this->value->removeEvent(time);
	if ( success == false )
	{
		throw RbException("Failed to remove time.");
	}
}

void MarkovTimesDistribution::addTime(double time)
{
	// remove the event
	bool success = this->value->addEvent(time);
	if ( success == false )
	{
		throw RbException("Failed to add time.");
	}
}

void MarkovTimesDistribution::simulateChildren()
{
	// check whether children need to be simulated
	const std::vector<DagNode*>& children = this->dag_node->getChildren();
	for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		if ( (*it)->isStochastic() == true )
		{
			AbstractEventsDistribution* dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
	        if ( dist != NULL )
	        {
	    		dist->resimulate();
	        }
		}
		else
		{
			(*it)->touch();
		}
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
