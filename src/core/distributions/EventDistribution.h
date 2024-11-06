#ifndef EventDistribution_H
#define EventDistribution_H

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements a generic mixture distribution between several possible values.
     *
     * This mixture can be considered as a multinomial distribution. We specify a vector of probabilities
     * and a vector of values. Then, a value drawn from this distribution takes each value corresponding to
     * its probability.
     * The values are already of the correct mixture type. You may want to apply a mixture allocation move
     * to change between the current value. The values themselves change automatically when the input parameters change.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    template <class valueType>
    class EventDistribution : public TypedDistribution< RbVector<valueType> > {
        
    public:
        // constructor(s)
        EventDistribution(TypedDistribution<std::int64_t> *ep, TypedDistribution<valueType> *vp);
        EventDistribution(const EventDistribution<valueType> &ep);

        // public member functions
        EventDistribution*                                  clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:

        // helper methods
        void                                                simulate();
        
        // private members
        TypedDistribution<std::int64_t>*                            event_prior;
        TypedDistribution<valueType>*                       value_prior;

    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>
#include <cstdint>

template <class valueType>
RevBayesCore::EventDistribution<valueType>::EventDistribution(TypedDistribution<std::int64_t> *ep, TypedDistribution<valueType> *vp) : TypedDistribution< RbVector<valueType> >( new RbVector<valueType>() ),
    event_prior( ep ),
    value_prior( vp )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& value_prior_pars = value_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    simulate();
}


template <class valueType>
RevBayesCore::EventDistribution<valueType>::EventDistribution( const EventDistribution<valueType> &d ) : TypedDistribution< RbVector<valueType> >(d),
    event_prior( d.event_prior->clone() ),
    value_prior( d.value_prior->clone() )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& value_prior_pars = value_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }    
    
}



template <class valueType>
RevBayesCore::EventDistribution<valueType>* RevBayesCore::EventDistribution<valueType>::clone( void ) const
{
    
    return new EventDistribution<valueType>( *this );
}



template <class valueType>
double RevBayesCore::EventDistribution<valueType>::computeLnProbability( void )
{
    size_t num_events = this->value->size();
    event_prior->setValue( new std::int64_t(num_events) );
    
    double ln_prob = event_prior->computeLnProbability();
    
    for (int i = 0; i < num_events; ++i)
    {
        
        // we also need to multiply with the probability of the value for this table
        value_prior->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( this->value->operator[](i) ) );
        ln_prob += value_prior->computeLnProbability();
        
    }
    
    return ln_prob;
}


template <class valueType>
void RevBayesCore::EventDistribution<valueType>::simulate()
{
    // clear the current value
    this->value->clear();
    
    // draw a number of events
    event_prior->redrawValue();
    std::int64_t num_events = event_prior->getValue();
    
    
    for (int i = 0; i < num_events; ++i)
    {
        
        // we also need to multiply with the probability of the value for this table
        value_prior->redrawValue();
        this->value->push_back( value_prior->getValue() );
        
    }
    
}


template <class valueType>
void RevBayesCore::EventDistribution<valueType>::redrawValue( void )
{
    
    simulate();
    
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::EventDistribution<valueType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    event_prior->swapParameter(oldP,newP);
    value_prior->swapParameter(oldP,newP);
    
}

#endif

