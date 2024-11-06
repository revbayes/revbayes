 /**
 * @file
 * This file contains the declaration and implementation
 * of the EmpiricalDistribution class. The distribution
 * is constructed from an input tree trace.
 * 
 *
 * @brief Declaration and implementation of the EmpiricalDistribution class
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0
 *
 */

#ifndef EmpiricalDistribution_H
#define EmpiricalDistribution_H

#include <cstdint>

#include "MemberObject.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDistribution.h"
#include "Trace.h"
#include "Tree.h"

namespace RevBayesCore {
    
    template<class valueType>
    class EmpiricalDistribution : public TypedDistribution<valueType>, public MemberObject<std::int64_t> {
        
    public:
		
		EmpiricalDistribution(Trace<valueType>* t, Trace<double>* d = NULL);
		
		virtual                                             ~EmpiricalDistribution(void);

		EmpiricalDistribution*                              clone(void) const;                                          //!< Create an independent clone
		double                                              computeLnProbability(void);                                 //!< Compute ln prob of current value
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;     //!< Map the member methods to internal function calls
        size_t                                              getBurnin(void) const;
        size_t                                              getCurrentIndex(void) const;
        size_t                                              getNumberOfSamples(void) const;
        void                                                redrawValue(void);                                          //!< Draw a new random value from distribution
		void                                                setCurrentIndex(size_t index);
        void                                                setValue(valueType *v, bool f=false);
		
	protected:

		void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
		
	private:
		
		size_t                                              current_index;
		Trace<double>*                                      density;
		Trace<valueType>*                                   trace;
		
	};

} 


template <class valueType>
RevBayesCore::EmpiricalDistribution<valueType>::EmpiricalDistribution(Trace<valueType>* t, Trace<double>* d) : TypedDistribution<valueType>( new valueType() ),
    density( d ),
    trace( t )
{
    // draw the first value
    redrawValue();
}


template <class valueType>
RevBayesCore::EmpiricalDistribution<valueType>::~EmpiricalDistribution()
{

}


template <class valueType>
RevBayesCore::EmpiricalDistribution<valueType>* RevBayesCore::EmpiricalDistribution<valueType>::clone( void ) const
{

    return new EmpiricalDistribution( *this );
}


template <class valueType>
void RevBayesCore::EmpiricalDistribution<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::int64_t &rv) const
{

    if ( n == "getIndex" )
    {
        rv = std::int64_t(current_index) + 1;
    }
    else
    {
        throw RbException("A empirical distribution does not have a member method called '" + n + "'.");
    }

}


template <class valueType>
size_t RevBayesCore::EmpiricalDistribution<valueType>::getBurnin( void ) const
{

    return trace->getBurnin();

}


template <class valueType>
size_t RevBayesCore::EmpiricalDistribution<valueType>::getCurrentIndex( void ) const
{

    return current_index;

}


template <class valueType>
size_t RevBayesCore::EmpiricalDistribution<valueType>::getNumberOfSamples( void ) const
{

    return trace->size();

}


template <class valueType>
void RevBayesCore::EmpiricalDistribution<valueType>::redrawValue( void )
{

    // draw a random tree
    RandomNumberGenerator* rng = GLOBAL_RNG;

    current_index = (size_t)( rng->uniform01() * trace->size(true) + trace->getBurnin() );

    delete this->value;
    this->value = new valueType( trace->objectAt(current_index) );;

}


template <class valueType>
void RevBayesCore::EmpiricalDistribution<valueType>::setCurrentIndex( size_t index )
{

    current_index = index;

    delete this->value;

    this->value = new valueType( trace->objectAt(current_index) );;

}


template <class valueType>
void RevBayesCore::EmpiricalDistribution<valueType>::setValue(valueType *v, bool f)
{

    bool found = false;
    for (size_t i = trace->getBurnin(); i < trace->size(); ++i)
    {
        if (trace->objectAt(i) == *v)
        {
            found = true;
            current_index = i;
            break;
        }
    }

    if (found == false)
    {
        RbException("The starting value is not in the empirical sample.");
    }

    TypedDistribution<Tree>::setValue(v, f);

}


template <class valueType>
double RevBayesCore::EmpiricalDistribution<valueType>::computeLnProbability( void )
{
    double lnProb = -log( trace->size(true) );

    if( density != NULL )
    {
        lnProb += density->objectAt(current_index);
    }

    return -log( trace->size(true) );

}


template <class valueType>
void RevBayesCore::EmpiricalDistribution<valueType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
}

#endif
