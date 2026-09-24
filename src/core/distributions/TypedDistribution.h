/**
 * @file
 * This file contains the declaration of the interface for all distributions.
 *
 * Distributions are either values inside DAG nodes, i.e. a constant node used as input for the DPP,
 * or be associated with a stochastic node.
 *
 * First, some distributions are requiring a distribution as a parameter, e.g. a generating distribution. Thus,
 * we need to implement distributions as objects storable in DAG nodes.
 *
 * Second, all stochastic nodes hold a distribution pointer. The value of the stochastic node is returned via
 * a call to get value in the distribution.
 *
 * Every distribution owns its value and hence this class is templated. Owning the value
 * has the advantage that calls to update can modify the value instead of creating a new object.
 * This is benefitial in functions generating large objects.
 *
 * @brief Declaration of distributions.
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-06-20 22:57:09 +0200 (Wed, 20 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-09-04, version 1.0
 *
 * $Id: Function.h 1643 2012-06-20 20:57:09Z hoehna $
 */


#ifndef TypedDistribution_H
#define TypedDistribution_H

#include "Distribution.h"
#include "Function.h"
#include "RbVector.h"

#include <iostream>

namespace RevBayesCore {
    
    template <class variableType>
    class StochasticNode;
    
    template<class variableType>
    class TypedDistribution : public Distribution {
        
    public:
        // constructors and destructor
        virtual                                        ~TypedDistribution(void);
               
        // public methods
        virtual const RevBayesCore::RbVector<variableType>&         getParameterValues(void) const;
        variableType&                                               getValue(void);                                                             //!< Get the current value (non-const)
        const variableType&                                         getValue(void) const;                                                       //!< Get the current value
        StochasticNode<variableType>*                               getStochasticNode(void);                                                    //!< Get the stochastic node holding this distribution
        void                                                        setOwnsValue(bool tf);                                                      //!< Set if we own the current value
        virtual void                                                setStochasticNode(StochasticNode<variableType> *n);                         //!< Set the stochastic node holding this distribution
        
        // virtual methods
        virtual void                                                setValue(variableType *v, bool f=false);                                    //!< Set the current value, e.g. attach an observation (clamp)
        virtual bool                                                allowsSA(void) { return false; }                                            //!< Checks if distribution is compatible with sampled ancestors
        
        // pure virtual public methods
        virtual TypedDistribution*                                  clone(void) const = 0;                                                      //!< Clone the distribution
        virtual double                                              computeLnProbability(void) = 0;                                             //!< Clone the ln probability density
        virtual void                                                redrawValue(SimulationCondition c);                                         //!< Draw a new random value from the distribution
        virtual void                                                redrawValue(void) = 0;                                                      //!< Draw a new random value from the distribution

    protected:
        TypedDistribution(variableType *v);
        TypedDistribution(const TypedDistribution &d);
        
        // overloaded operators
        TypedDistribution&                                          operator=(const TypedDistribution &d);

        virtual void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;        //!< Exchange the parameter

        
        // inheritable attributes
        StochasticNode<variableType>*                               dag_node;                                                                   //!< The stochastic node holding this distribution. This is needed for delegated calls to the DAG, such as getAffected(), ...
        bool                                                        owns_value;
        variableType*                                               value;
        
    };
    
    // Global functions using the class
    template <class variableType>
    std::ostream&                                                   operator<<(std::ostream& o, const TypedDistribution<variableType>& x);                                //!< Overloaded output operator
    
}


#include "Cloner.h"
#include "Cloneable.h"
#include "IsDerivedFrom.h"
#include "RbException.h"


template <class variableType>
RevBayesCore::TypedDistribution<variableType>::TypedDistribution(variableType *v) : Distribution(), 
    dag_node( NULL ),
    value( v ),
    owns_value( true )
{
    
}

template <class variableType>
RevBayesCore::TypedDistribution<variableType>::TypedDistribution(const TypedDistribution &d) : Distribution(d), 
    dag_node( NULL ),
    value( d.owns_value ? Cloner<variableType, IsDerivedFrom<variableType, Cloneable>::Is >::createClone( *d.value ) : d.value ),
    owns_value( d.owns_value )
{
    
}

template <class variableType>
RevBayesCore::TypedDistribution<variableType>::~TypedDistribution( void )
{
    
    if ( owns_value == true )
    {
        delete value;
    }
    
}



template <class variableType>
RevBayesCore::TypedDistribution<variableType>& RevBayesCore::TypedDistribution<variableType>::operator=(const TypedDistribution &d)
{
    
    if ( this != &d ) 
    {
        // call base class
        Distribution::operator=( d );

        // clean up if necessary
        if ( owns_value == true )
        {
            delete value;
        }
        
        // copy member variables
        owns_value = d.owns_value;
        
        // make my own copy of the value (we rely on proper implementation of assignment operators)
        value = ( owns_value ? Cloner<variableType, IsDerivedFrom<variableType, Cloneable>::Is >::createClone( *d.value ) : d.value );
    }
    
    return *this;
}


template <class variableType>
const RevBayesCore::RbVector<variableType>& RevBayesCore::TypedDistribution<variableType>::getParameterValues(void) const
{
    throw RbException("Cannot access mixture values of non-mixture distribution.");
}

template <class variableType>
variableType& RevBayesCore::TypedDistribution<variableType>::getValue(void)
{
    
    return *value;
}

template <class variableType>
const variableType& RevBayesCore::TypedDistribution<variableType>::getValue(void) const
{
    
    return *value;
}


template <class variableType>
RevBayesCore::StochasticNode<variableType>* RevBayesCore::TypedDistribution<variableType>::getStochasticNode( void )
{
    
    return dag_node;
}

template <class variableType>
void RevBayesCore::TypedDistribution<variableType>::redrawValue(SimulationCondition c)
{
    
    // by default we delegate to the method without arguments
    // this allows us to overload this function but only if the argument is necessary
    redrawValue();
    
}

template <class variableType>
void RevBayesCore::TypedDistribution<variableType>::setOwnsValue(bool tf)
{
    
    // create our own copy of the value if we didn't own it before
    if ( value != NULL && owns_value == false && tf == true )
    {
        value = Cloner<variableType, IsDerivedFrom<variableType, Cloneable>::Is >::createClone( *value );
    }

    // delete the value if we are not supposed to own it anymore
    if ( value != NULL && owns_value == true && tf == false )
    {
        delete value;
    }
    
    owns_value = tf;
}

template <class variableType>
void RevBayesCore::TypedDistribution<variableType>::setStochasticNode( StochasticNode<variableType> *n )
{
    
    dag_node = n;
}

template <class variableType>
void RevBayesCore::TypedDistribution<variableType>::setValue( variableType *v, bool force )
{
    
    // free memory
    if ( value != v && owns_value == true )
    {
        delete value;
    }
    
    if ( force == false && owns_value == false )
    {
        (*value) = (*v);
        delete v;
    }
    else
    {
        value = v;
    }
        
}


template <class variableType>
std::ostream& RevBayesCore::operator<<(std::ostream& o, const TypedDistribution<variableType>& f)
{
    
    o << "Distribution()";
    
    return o;
}

#endif
