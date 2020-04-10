#ifndef MarkovEventsDistribution_H
#define MarkovEventsDistribution_H

#include <vector>

#include "OrderedEvents.h"
#include "RbVector.h"
#include "OrderedEventTimes.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    template <class valueType>
    class MarkovEventsDistribution : public TypedDistribution< OrderedEvents<valueType> > {
        
    public:
        // constructor(s)
        MarkovEventsDistribution(const TypedDagNode< OrderedEventTimes > *oet, TypedDistribution<valueType> *g);
        MarkovEventsDistribution(const MarkovEventsDistribution &d);
        
        // public member functions
        MarkovEventsDistribution*                           clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // helper methods
        OrderedEvents<valueType>*                           simulate(void);
              
        // private members
        const TypedDagNode< OrderedEventTimes > *           event_times;
        TypedDistribution<valueType>*                       base_distribution;

    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"


/** MarkovEventsDistribution Constructor
 * @param g the base distribution
 * @param cp concentration parameter
 * @param n number of elements
 */
template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>::MarkovEventsDistribution(const TypedDagNode< OrderedEventTimes > *oet, TypedDistribution<valueType> *g) :
                                                                            TypedDistribution< OrderedEvents<valueType> >( new OrderedEvents<valueType>() ),
    event_times( oet ),
    base_distribution( g )
{
    this->addParameter( event_times );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }

    this->value = simulate();
}


template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>::MarkovEventsDistribution(const MarkovEventsDistribution &d) : TypedDistribution< OrderedEvents<valueType> >( d ),
	event_times( d.event_times ),
	base_distribution( d.base_distribution )
{
    this->addParameter( event_times );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    this->value = simulate();
}



template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>* RevBayesCore::MarkovEventsDistribution<valueType>::clone( void ) const
{
    return new MarkovEventsDistribution<valueType>( *this );
}

template <class valueType>
double RevBayesCore::MarkovEventsDistribution<valueType>::computeLnProbability( void )
{
	double ln_prob = 0.0;

	const typename std::map<double, valueType> &the_values = this->value->getEvents();
	for( typename std::map<double, valueType>::const_iterator it = the_values.begin(); it != the_values.end(); ++it)
	{
		// set the value for the base distribution
		base_distribution->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( it->second ) );

		// compute the probability
		ln_prob += base_distribution->computeLnProbability();
	}

    return ln_prob;
}



template <class valueType>
RevBayesCore::OrderedEvents<valueType>* RevBayesCore::MarkovEventsDistribution<valueType>::simulate()
{
	// get the number of events from the times
	const std::set<double>& the_times  = event_times->getValue().getEventTimes();
	size_t                  num_events = event_times->getValue().size();

	// create the new container
	OrderedEvents<valueType>* sim = new OrderedEvents<valueType>();

	// for each event, simulate a value
	for(std::set<double>::const_iterator it = the_times.begin(); it != the_times.end(); ++it)
	{
		// have the base distribution simulate a new value
		base_distribution->redrawValue();

		// copy the simulated value
		valueType new_event = base_distribution->getValue();

		// add the new event
		sim->addEvent(*it, new_event);
	}

	// return the simulation
	return sim;
}


template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::redrawValue( void )
{
    delete this->value;
    this->value = simulate();
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == event_times)
    {
        event_times = static_cast<const TypedDagNode< OrderedEventTimes > *>( newP );
    }
    else
    {
        base_distribution->swapParameter(oldP,newP);
    }
}




#endif

