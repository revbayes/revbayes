/**
 * @file
 * This file contains the declaration of the interface for all functions.
 *
 * Functions either be values inside DAG nodes, i.e. the random variable of some distribution,
 * or be associated with a deterministic node.
 *
 * First, some distributions are requiring a function as a parameter, e.g. a rate function. Thus,
 * we need to implement functions as objects storable in DAG nodes.
 *
 * Second, all deterministic nodes hold a function. The value of the deterministic node is returned via
 * a call to get value in the function (and potentially delegate to update the value if it is dirty).
 *
 * Every functions owns it's value and hence this class is templated. Owning the value
 * has the advantage that calls to update can modify the value instead of creating a new object.
 * This is benefitial in functions generating large objects.
 *
 * @brief Declaration of functions.
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-06-20 22:57:09 +0200 (Wed, 20 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-09-02, version 1.0
 *
 * $Id: Function.h 1643 2012-06-20 20:57:09Z hoehna $
 */

#ifndef TypedFunction_H
#define TypedFunction_H

#include "Function.h"

#include <iostream>

namespace RevBayesCore {
    
    template <class valueType>
    class DeterministicNode;
    
    template<class valueType>
    class TypedFunction : public Function {
        
    public:
        // constructors and destructor
        virtual                            ~TypedFunction(void);
        
        // public methods
        virtual valueType&                  getValue(void);                                                             //!< Get a value reference
        virtual const valueType&            getValue(void) const;                                                       //!< Get value reference (const)
        void                                setDeterministicNode(DeterministicNode<valueType> *n);                      //!< Set the stochastic node holding this distribution

        // pure virtual public methors
        virtual TypedFunction*              clone(void) const = 0;                                                      //!< Clone the function
        virtual void                        touch(const RevBayesCore::DagNode *toucher );                               //!< Touch the function (set dirty flag)
        virtual void                        update(void) = 0;                                                           //!< Update the value of the function
    
    protected:
        TypedFunction(valueType *v);
        TypedFunction(const TypedFunction &f);
       
        // overloaded operators
        TypedFunction&                      operator=(const TypedFunction &d); 
        
        virtual void                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;        //!< Exchange the parameter

        // members 
        DeterministicNode<valueType>*       dag_node;                                                                    //!< The deterministic node holding this function. This is needed for delegated calls to the DAG, such as getAffected(), addTouchedElementIndex()...
        valueType*                          value;
    
    };
    
    // Global functions using the class
    template <class valueType>
    std::ostream&                           operator<<(std::ostream& o, const TypedFunction<valueType>& x);                                //!< Overloaded output operator
    
}

#include "Cloneable.h"
#include "Cloner.h"
#include "IsDerivedFrom.h"

template <class valueType>
RevBayesCore::TypedFunction<valueType>::TypedFunction(valueType *v) : Function(),
    dag_node( NULL ),
    value( v )
{
    
}

template <class valueType>
RevBayesCore::TypedFunction<valueType>::TypedFunction(const TypedFunction &f) : Function(f), 
    dag_node( NULL ),
    value( NULL )
{
    
    if ( f.value != NULL )
    {
        value = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( *f.value );
    }
    
}

template <class valueType>
RevBayesCore::TypedFunction<valueType>::~TypedFunction( void )
{
    
    delete value;
}


template <class valueType>
RevBayesCore::TypedFunction<valueType>& RevBayesCore::TypedFunction<valueType>::operator=(const TypedFunction &f)
{
    
    if ( this != &f ) 
    {
        // call base class
        Function::operator=( f );
        
        // make my own copy of the value (we rely on proper implementation of assignment operators)
        delete value;
        value = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( *f.value );
        
    }
    
    return *this;
}


template <class valueType>
const valueType& RevBayesCore::TypedFunction<valueType>::getValue(void) const 
{
        
    return *value;
}

template <class valueType>
valueType& RevBayesCore::TypedFunction<valueType>::getValue(void) 
{
        
    return *value;
}


template <class valueType>
void RevBayesCore::TypedFunction<valueType>::setDeterministicNode(DeterministicNode<valueType> *n) 
{
    
    dag_node = n;
    
}


template <class valueType>
void RevBayesCore::TypedFunction<valueType>::touch( const RevBayesCore::DagNode* /*toucher*/ )
{
    // nothing to do here in the base class

}


template <class valueType>
std::ostream& RevBayesCore::operator<<(std::ostream& o, const TypedFunction<valueType>& f) 
{
    
    o << "f(x) = " << f.getValue();
    
    return o;
}


#endif
