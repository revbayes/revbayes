#ifndef OrderedEventScaleProposal_H
#define OrderedEventScaleProposal_H

#include <set>
#include <string>

#include "OrderedEvents.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    template <class valueType>
    class OrderedEventScaleProposal : public Proposal {
        
    public:
        OrderedEventScaleProposal( StochasticNode< OrderedEvents<valueType> > *n, double l=0.5 );                                                                    //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        OrderedEventScaleProposal*          clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                              getProposalTuningParameter(void) const;
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                setProposalTuningParameter(double tp);
        void                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:

        // parameters
        StochasticNode< OrderedEvents<valueType> >*        variable;                                                                           //!< The variable the Proposal is working on
        double                                             lambda;
        valueType                                          old_value;
        double                                             event_time;
        bool                                               abort;
        
    };
    
}


#include "MixtureDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "MarkovEventsDistribution.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template <class valueType>
RevBayesCore::OrderedEventScaleProposal<valueType>::OrderedEventScaleProposal( StochasticNode< OrderedEvents<valueType> > *n, double l ) : Proposal(),
    variable( n ),
    lambda( l ),
	old_value( 0.0 ),
	event_time( 0.0 ),
	abort( false )

{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class valueType>
RevBayesCore::OrderedEventScaleProposal<valueType>* RevBayesCore::OrderedEventScaleProposal<valueType>::clone( void ) const
{
    
    return new OrderedEventScaleProposal<valueType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class valueType>
const std::string& RevBayesCore::OrderedEventScaleProposal<valueType>::getProposalName( void ) const
{

//	valueType::

    static std::string name = "OrderedEventScale";
    
    return name;
}


template <class valueType>
double RevBayesCore::OrderedEventScaleProposal<valueType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * The reversible jump proposal switches the current "dimension".
 *
 * \return The hastings ratio.
 */
template <class valueType>
double RevBayesCore::OrderedEventScaleProposal<valueType>::doProposal( void )
{
	// clear abort flag
	abort = false;

    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the distribution
    TypedEventsDistribution<valueType>* dist = dynamic_cast<TypedEventsDistribution<valueType> *>( &variable->getDistribution() );
    if ( dist == NULL )
    {
    	throw RbException("Tried to use OrderedEventScaleProposal on an invalid type.");
    }

    // get the random event
    // abort if there are no events
    double this_event_time;
    valueType this_event_value;
    abort = !dist->getRandomEvent(this_event_time, this_event_value);
    if ( abort == true )
    {
    	return RbConstants::Double::neginf;
    }

    // get the value
    event_time = this_event_time;
    old_value  = this_event_value;

    // perform a Scale proposal
	double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
	double new_value = old_value * scaling_factor;

	// update the value in the variable
	abort = !variable->getValue().changeEvent(event_time, new_value);
	if ( abort == true )
	{
		return RbConstants::Double::neginf;
	}

	// touch the variable
	variable->touch(true);

	return std::log(scaling_factor);
}


/**
 *
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    // nothing to print
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::undoProposal( void )
{
	if ( abort == false )
	{
	    variable->getValue().changeEvent(event_time, old_value);
	    variable->touch(true);
	}
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<OrderedEvents<valueType> >* >(newN) ;
    
}


template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::setProposalTuningParameter(double tp)
{
	lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
template <class valueType>
void RevBayesCore::OrderedEventScaleProposal<valueType>::tune( double rate )
{
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
    	lambda /= (2.0 - rate/p);
    }
    
    lambda = fmin(10000, lambda);
}



#endif

