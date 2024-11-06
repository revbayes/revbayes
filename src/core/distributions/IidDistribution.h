#ifndef IidDistribution_H
#define IidDistribution_H

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
    class IidDistribution : public TypedDistribution< RbVector<valueType> > {
        
    public:
        // constructor(s)
        IidDistribution(std::int64_t n, TypedDistribution<valueType> *vp);
        IidDistribution(const IidDistribution<valueType> &d);

        // public member functions
        IidDistribution*                                    clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        
        // helper methods
        void                                                simulate();
        
        // private members
        std::int64_t                                                n_samples;
        TypedDistribution<valueType>*                       value_prior;
        
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>
#include <cstdint>

template <class valueType>
RevBayesCore::IidDistribution<valueType>::IidDistribution(std::int64_t n, TypedDistribution<valueType> *vp) : TypedDistribution< RbVector<valueType> >( new RbVector<valueType>() ),
    n_samples( n ),
    value_prior( vp )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& value_prior_pars = value_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    simulate();
}


template <class valueType>
RevBayesCore::IidDistribution<valueType>::IidDistribution( const IidDistribution<valueType> &d ) : TypedDistribution< RbVector<valueType> >(d),
    n_samples( d.n_samples ),
    value_prior( d.value_prior->clone() )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& value_prior_pars = value_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }

    
}



template <class valueType>
RevBayesCore::IidDistribution<valueType>* RevBayesCore::IidDistribution<valueType>::clone( void ) const
{
    
    return new IidDistribution<valueType>( *this );
}



template <class valueType>
double RevBayesCore::IidDistribution<valueType>::computeLnProbability( void )
{
    
    double ln_prob = 0.0;
    
    for (int i = 0; i < n_samples; ++i)
    {

        // we also need to multiply with the probability of the value for this table
        value_prior->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( this->value->operator[](i) ) );
        ln_prob += value_prior->computeLnProbability();
        
    }

    return ln_prob;
}


template <class valueType>
void RevBayesCore::IidDistribution<valueType>::simulate()
{
    // clear the current value
    this->value->clear();

    for (int i = 0; i < n_samples; ++i)
    {

        // we also need to multiply with the probability of the value for this table
        value_prior->redrawValue();
        this->value->push_back( value_prior->getValue() );

    }
    
}


template <class valueType>
void RevBayesCore::IidDistribution<valueType>::redrawValue( void )
{
    
    simulate();
    
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::IidDistribution<valueType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    value_prior->swapParameter(oldP,newP);
    
}

#endif

