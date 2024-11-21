#ifndef ReversibleJumpMixtureConstantDistribution_H
#define ReversibleJumpMixtureConstantDistribution_H

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements a generic mixture distribution between a constant value and a distribution.
     *
     * The mixture consists of two components:
     *   a) a constant value
     *   b) a distribution
     * with probability p the current value is from a) and with probability 1-p it is from b).
     * Here we need reversible jump MCMC to mix/jump between the two possible states.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    template <class mixtureType>
    class ReversibleJumpMixtureConstantDistribution : public TypedDistribution<mixtureType>, public MemberObject<long> {
        
    public:
        // constructor(s)
        ReversibleJumpMixtureConstantDistribution(const TypedDagNode< mixtureType > *cv, TypedDistribution< mixtureType > *dv, const TypedDagNode<double> *p);
        ReversibleJumpMixtureConstantDistribution(const ReversibleJumpMixtureConstantDistribution<mixtureType> &d);

        virtual                                            ~ReversibleJumpMixtureConstantDistribution();
        
        ReversibleJumpMixtureConstantDistribution<mixtureType>&     operator=(const ReversibleJumpMixtureConstantDistribution<mixtureType> &d);

        
        // public member functions
        ReversibleJumpMixtureConstantDistribution*                  clone(void) const;                                                                              //!< Create an independent clone
        double                                                      computeLnProbability(void);
        void                                                        executeMethod(const std::string &n, const std::vector<const DagNode*> &args, long &rv) const;    //!< Map the member methods to internal function calls
        const TypedDistribution<mixtureType>&                       getBaseDistribution(void) const;
        TypedDistribution<mixtureType>&                             getBaseDistribution(void);
        const mixtureType&                                          getConstantValue(void) const;
        size_t                                                      getCurrentIndex(void) const;
        size_t                                                      getNumberOfCategories(void) const;
        void                                                        redrawValue(void);
        void                                                        redrawValueByIndex(int i);
        void                                                        setCurrentIndex(size_t i);
        void                                                        setValue(mixtureType *v, bool f=false);
        
        // special handling of state changes
        void                                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                          //!< get affected nodes
        void                                                        keepSpecialization(const DagNode* affecter);
        void                                                        restoreSpecialization(const DagNode *restorer);
        void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);

    protected:
        // Parameter management functions
        void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
        
    private:
        // helper methods
        mixtureType*                                                simulate();
        
        // private members
        const TypedDagNode< mixtureType >*                          const_value;
        TypedDistribution<mixtureType>*                                base_distribution;
        const TypedDagNode< double >*                               probability;
        
        size_t                                                      index;
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>

template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::ReversibleJumpMixtureConstantDistribution(const TypedDagNode< mixtureType > *cv, TypedDistribution< mixtureType > *dv, const TypedDagNode< double > *p) : TypedDistribution<mixtureType>( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( cv->getValue() ) ),
    const_value( cv ),
    base_distribution( dv ),
    probability( p ),
    index( 0 )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( const_value );
    this->addParameter( probability );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    delete this->value;
    this->value = simulate();
}


template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::ReversibleJumpMixtureConstantDistribution(const ReversibleJumpMixtureConstantDistribution<mixtureType> &d) : TypedDistribution<mixtureType>( d ),
    const_value( d.const_value ),
    base_distribution( d.base_distribution->clone() ),
    probability( d.probability ),
    index( d.index )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( const_value );
    this->addParameter( probability );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
}



template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>& RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::operator=(const RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType> &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<mixtureType>::operator=( d );
        
        delete base_distribution;
        
        const_value         = d.const_value;
        base_distribution   = d.base_distribution->clone();
        probability         = d.probability;
        index               = d.index;
        
    }
    
    return *this;
}

template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::~ReversibleJumpMixtureConstantDistribution( void )
{

    delete base_distribution;
    
}


template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>* RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::clone( void ) const
{
    
    return new ReversibleJumpMixtureConstantDistribution<mixtureType>( *this );
}



template <class mixtureType>
double RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::computeLnProbability( void )
{
    
    double ln_prob;
    if ( index == 0 )
    {
        if ( *this->value != const_value->getValue() )
        {
            ln_prob = RbConstants::Double::neginf;
        }
        else
        {
            ln_prob = log( probability->getValue() );
        }
        
    }
    else
    {
        
        ln_prob = log( 1.0 - probability->getValue() );
        base_distribution->setValue( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone(*this->value) );
        ln_prob += base_distribution->computeLnProbability();
        
    }
    
    return ln_prob;
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, long &rv) const
{
    
    if ( n == "index" )
    {
        rv = long(index);
    }
    else
    {
        throw RbException() << "A mixture distribution does not have a member method called '" << n << "'.";
    }
    
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter)
{
    // only do this when the toucher was our constant value and this value was supposed to be equal to the constant value
    if ( affecter == const_value && index == 0 && this->dag_node != NULL )
    {
        this->dag_node->initiateGetAffectedNodes( affected );
    }
    
}


template <class mixtureType>
const RevBayesCore::TypedDistribution<mixtureType>& RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getBaseDistribution( void ) const
{
    
    return *base_distribution;
}


template <class mixtureType>
RevBayesCore::TypedDistribution<mixtureType>& RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getBaseDistribution( void )
{
    
    return *base_distribution;
}


template <class mixtureType>
const mixtureType& RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getConstantValue( void ) const
{
    
    return const_value->getValue();
}



template <class mixtureType>
size_t RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getCurrentIndex( void ) const
{
    
    return index;
}


template <class mixtureType>
size_t RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::getNumberOfCategories( void ) const
{
    
    return 2;
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::keepSpecialization( const DagNode* affecter )
{
    // only do this when the toucher was our constant value and this value was supposed to be equal to the constant value
    if ( affecter == const_value && index == 0 && this->dag_node != NULL )
    {
        this->dag_node->keepAffected();
    }
    
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::redrawValue( void )
{
    
    delete this->value;
    this->value = simulate();
    
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::redrawValueByIndex( int i )
{
    delete this->value;
    
    if (i == 0)
    {
        index = 0;
        this->value = Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( const_value->getValue() );
    }
    else
    {
        index = 1;
        base_distribution->redrawValue();
        this->value = Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( base_distribution->getValue() );
    }
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::restoreSpecialization( const DagNode *restorer )
{
    // only do this when the toucher was our constant value and this value was supposed to be equal to the constant value
    if ( restorer == const_value && index == 0 )
    {
        if constexpr (std::is_base_of_v<Cloneable,mixtureType>)
        {
            delete this->value;
            this->value = const_value->getValue().clone();
        }
        else
        {
            *this->value = const_value->getValue();
        }

        if ( this->dag_node != NULL )
        {
            this->dag_node->restoreAffected();
        }
    }
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::setCurrentIndex(size_t i)
{
    // we complete rely on the move to also change the value accordingly!!!
    index = i;
}


template <class mixtureType>
mixtureType* RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::simulate()
{
    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    
    if ( u < probability->getValue() )
    {
        index = 0;
        return Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( const_value->getValue() );
    }
    else
    {
        index = 1;
        base_distribution->redrawValue();
        return Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( base_distribution->getValue() );
    }
    
}


/** Swap a parameter of the distribution */
template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == const_value)
    {
        const_value = static_cast<const TypedDagNode< mixtureType >* >( newP );
    }
    else if (oldP == probability)
    {
        probability = static_cast<const TypedDagNode< double >* >( newP );
    }
    else
    {
        base_distribution->swapParameter(oldP,newP);
    }
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::setValue(mixtureType *v, bool force)
{
    
    if(v != this->value) delete this->value;
    
    if ( force == false )
    {
        this->value = v;
    }
    else
    {
        
        if ( *v == const_value->getValue() )
        {
            index = 0;
            this->value = v;
        }
        else
        {
            index = 1;
            base_distribution->setValue( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( *v ) );
            this->value = v;
        }
        
    }
    
}


template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureConstantDistribution<mixtureType>::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    // only do this when the toucher was our constant value and this value was supposed to be equal to the constant value
    if ( toucher == const_value && index == 0 )
    {
        if constexpr (std::is_base_of_v<Cloneable,mixtureType>)
        {
            delete this->value;
            this->value = const_value->getValue().clone();
        }
        else
        {
            *this->value = const_value->getValue();
        }
        
        if ( this->dag_node != NULL )
        {
            this->dag_node->touchAffected();
        }
    }
}

#endif
