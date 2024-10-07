#ifndef ProbabilityDensityFunction_H
#define ProbabilityDensityFunction_H

#include "ContinuousFunction.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * This class implements a generic probability density function of a distribution.
     *
     *
     * The distribution is given as a parameter so this file is the wrapper to call the pdf of the distribution.
     * Hence, this function can be used inside deterministic nodes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-12-1, version 1.0
     */
    template <class valueType>
    class ProbabilityDensityFunction : public ContinuousFunction {
        
    public:
        ProbabilityDensityFunction(const TypedDagNode<valueType> *x, TypedDistribution<valueType> *d, bool l=true);
        ProbabilityDensityFunction(const ProbabilityDensityFunction &d);
        virtual                            ~ProbabilityDensityFunction(void);
        
        ProbabilityDensityFunction&         operator=(const ProbabilityDensityFunction &d);

        
        ProbabilityDensityFunction*         clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        
        const TypedDagNode<valueType>*      x;
        TypedDistribution<valueType>*       dist;
        bool                                useLog;
        
    };
}

template <class valueType>
RevBayesCore::ProbabilityDensityFunction<valueType>::ProbabilityDensityFunction(const TypedDagNode<valueType> *z, TypedDistribution<valueType>* d, bool l) : ContinuousFunction( new double(0.0) ),
    x( z ),
    dist( d ),
    useLog( l )
{

    addParameter( x );
    
}


template <class valueType>
RevBayesCore::ProbabilityDensityFunction<valueType>::ProbabilityDensityFunction(const ProbabilityDensityFunction &d) : ContinuousFunction( d ),
    x( d.x ),
    dist( d.dist->clone() ),
    useLog( d.useLog )
{
    
    addParameter( x );
    
}

template <class valueType>
RevBayesCore::ProbabilityDensityFunction<valueType>::~ProbabilityDensityFunction(void)
{
    delete dist;
}


template <class valueType>
RevBayesCore::ProbabilityDensityFunction<valueType>& RevBayesCore::ProbabilityDensityFunction<valueType>::operator=(const ProbabilityDensityFunction &d)
{
    
   // check for self-assignment
    if ( this != d )
    {
        ContinuousFunction::operator=( d );
        delete dist;
        
        x       = d.x;
        dist    = d.dist->clone();
        useLog  = d.useLog;

    }
    
    return *this;
}


template <class valueType>
RevBayesCore::ProbabilityDensityFunction<valueType>* RevBayesCore::ProbabilityDensityFunction<valueType>::clone( void ) const
{
    return new ProbabilityDensityFunction<valueType>(*this);
}


template <class valueType>
void RevBayesCore::ProbabilityDensityFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<valueType>* >( newP );
    }
    else
    {
        dist->swapParameter(oldP, newP);
    }
    
}

template <class valueType>
void RevBayesCore::ProbabilityDensityFunction<valueType>::update( void )
{
     if constexpr(std::is_base_of_v<Cloneable, valueType>) {
        // For abstract classes that cannot be instantiated directly
        // (e.g., AbstractHomologousDiscreteCharacterData), we clone the
        // object and pass a pointer to the cloned value
        dist->setValue(const_cast<valueType*>(&x->clone()->getValue()));
    } else {
        // For other types, we create a new instance
        dist->setValue(new valueType(x->getValue()));
    }
    
    if ( useLog == true )
    {
        *value = dist->computeLnProbability();
    }
    else
    {
        *value = exp( dist->computeLnProbability() );
    }
    
}

#endif
