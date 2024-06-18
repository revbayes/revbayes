#ifndef ShiftEventsFunction_H
#define ShiftEventsFunction_H

#include "OrderedEvents.h"
#include "RbVector.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class ShiftEventsFunction : public TypedFunction< OrderedEvents< valueType > > {
        
    public:
        ShiftEventsFunction(const TypedDagNode< valueType >* iv, const TypedDagNode< OrderedEvents<valueType> > *se);
        virtual                                            ~ShiftEventsFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        ShiftEventsFunction*                                clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< valueType >*                    initial_values;
        const TypedDagNode< OrderedEvents<valueType> >*     shift_events;
        
    };
    
}




template <class valueType>
RevBayesCore::ShiftEventsFunction<valueType>::ShiftEventsFunction(const TypedDagNode< valueType >* iv, const TypedDagNode<OrderedEvents<valueType>> *se) : TypedFunction< OrderedEvents< valueType > >( new OrderedEvents<valueType>() ),
	initial_values( iv ),
	shift_events( se )
{
    
    // add the parameter as a parent
    this->addParameter( initial_values );
    this->addParameter( shift_events );
    
    update();
}


template <class valueType>
RevBayesCore::ShiftEventsFunction<valueType>::~ShiftEventsFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::ShiftEventsFunction<valueType>* RevBayesCore::ShiftEventsFunction<valueType>::clone( void ) const
{
    return new ShiftEventsFunction<valueType>( *this );
}



template <class valueType>
void RevBayesCore::ShiftEventsFunction<valueType>::update( void )
{
	// empty the current state
	delete this->value;

	// create the new object
	OrderedEvents<valueType>* new_value = new OrderedEvents<valueType>();

	// add the initial value
	valueType current_value = initial_values->getValue();
	new_value->addEvent(0.0, current_value );

	// add the shifts
	const std::map<double, valueType>& events = shift_events->getValue().getEvents();
	for(typename std::map<double, valueType>::const_iterator it = events.begin(); it != events.end(); ++it)
	{
		// multiply the current event
		current_value = current_value * it->second;

		// add the event
		new_value->addEvent(it->first, current_value);
	}

	// set the current state
	this->value = new_value;
}



template <class valueType>
void RevBayesCore::ShiftEventsFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == initial_values )
    {
    	initial_values = static_cast<const TypedDagNode<valueType>* >( newP );
        this->update();
    }
	else if ( oldP == shift_events )
    {
		shift_events = static_cast<const TypedDagNode<OrderedEvents<valueType>>* >( newP );
        this->update();
    }
    
}

#endif
