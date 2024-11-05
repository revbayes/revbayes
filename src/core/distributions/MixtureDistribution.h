#ifndef MixtureDistribution_H
#define MixtureDistribution_H

#include "MemberObject.h"
#include "RbVector.h"
#include "Simplex.h"
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
    template <class mixtureType>
    class MixtureDistribution : public TypedDistribution<mixtureType>, public MemberObject<std::int64_t> {
        
    public:
        // constructor(s)
        MixtureDistribution(const TypedDagNode< RbVector<mixtureType> > *v, const TypedDagNode< Simplex > *p);
        
        // public member functions
        MixtureDistribution*                                clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;     //!< Map the member methods to internal function calls
        const RevBayesCore::RbVector<mixtureType>&          getParameterValues(void) const;
        size_t                                              getCurrentIndex(void) const;
        std::vector<double>                                 getMixtureProbabilities(void) const;
        size_t                                              getNumberOfMixtureElements(void) const;                                                        //!< Get the number of elements for this value
        void                                                redrawValue(void);
        void                                                setCurrentIndex(size_t i);
        void                                                setValue(mixtureType *v, bool f=false);
        
        // special handling of state changes
        void                                                getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                          //!< get affected nodes
        void                                                keepSpecialization(const DagNode* affecter);
        void                                                restoreSpecialization(const DagNode *restorer);
        void                                                touchSpecialization(const DagNode *toucher, bool touchAll);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        // helper methods
        const mixtureType&                                  simulate();
        
        // private members
        const TypedDagNode< RbVector<mixtureType> >*        parameter_values;
        const TypedDagNode< Simplex >*                      probabilities;
        
        size_t                                              index;
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>

template <class mixtureType>
RevBayesCore::MixtureDistribution<mixtureType>::MixtureDistribution(const TypedDagNode< RbVector<mixtureType> > *v, const TypedDagNode< Simplex > *p) : TypedDistribution<mixtureType>( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( v->getValue()[0] ) ),
    parameter_values( v ),
    probabilities( p ),
    index( 0 )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( parameter_values );
    this->addParameter( probabilities );
    
    redrawValue();
}



template <class mixtureType>
RevBayesCore::MixtureDistribution<mixtureType>* RevBayesCore::MixtureDistribution<mixtureType>::clone( void ) const
{

    return new MixtureDistribution<mixtureType>( *this );
}



template <class mixtureType>
double RevBayesCore::MixtureDistribution<mixtureType>::computeLnProbability( void )
{
    
    return log(probabilities->getValue()[index]);
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::int64_t &rv) const
{
    
    if ( n == "getAllocationIndex" )
    {
        rv = long(index) + 1;
    }
    else
    {
        throw RbException("A mixture distribution does not have a member method called '" + n + "'.");
    }
    
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter)
{
    // only delegate when the toucher was our parameters
    if ( affecter == parameter_values && this->dag_node != NULL )
    {
        this->dag_node->initiateGetAffectedNodes( affected );
    }
    
}


template <class mixtureType>
size_t RevBayesCore::MixtureDistribution<mixtureType>::getCurrentIndex( void ) const
{

    return index;
}


template <class mixtureType>
std::vector<double> RevBayesCore::MixtureDistribution<mixtureType>::getMixtureProbabilities( void ) const
{

    return probabilities->getValue();
}


template <class mixtureType>
size_t RevBayesCore::MixtureDistribution<mixtureType>::getNumberOfMixtureElements( void ) const
{

    return probabilities->getValue().size();
}


template <class mixtureType>
const RevBayesCore::RbVector<mixtureType>& RevBayesCore::MixtureDistribution<mixtureType>::getParameterValues( void ) const
{
    return parameter_values->getValue();
}

template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::keepSpecialization( const DagNode* affecter )
{
    // only do this when the toucher was our parameters
    if ( affecter == parameter_values && this->dag_node != NULL )
    {
        this->dag_node->keepAffected();
    }
    
}


template <class mixtureType>
const mixtureType& RevBayesCore::MixtureDistribution<mixtureType>::simulate()
{
    
    const std::vector<double> &probs = probabilities->getValue();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    index = 0;
    while ( u > probs[index] )
    {
        u -= probs[index];
        ++index;
    }
    
    return parameter_values->getValue()[index];
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::redrawValue( void )
{
    if constexpr(std::is_base_of_v<Cloneable,mixtureType>)
    {
	delete this->value;
	this->value = simulate().clone();
    }
    else
	(*this->value) = simulate();

}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::setCurrentIndex(size_t i)
{
    index = i;

    const mixtureType &tmp = parameter_values->getValue()[i];

    if constexpr(std::is_base_of_v<Cloneable, mixtureType>)
    {
	delete this->value;
	this->value = tmp.clone();
    }
    else
	(*this->value) = tmp;
}


/** Swap a parameter of the distribution */
template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == parameter_values)
    {
        parameter_values = static_cast<const TypedDagNode< RbVector<mixtureType> >* >( newP );
    }
    else if (oldP == probabilities)
    {
        probabilities = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::restoreSpecialization( const DagNode *restorer )
{
    
    // only do this when the toucher was our parameters
    if ( restorer == parameter_values )
    {
        const mixtureType &tmp = parameter_values->getValue()[index];
        if constexpr(std::is_base_of_v<Cloneable, mixtureType>)
        {
            delete this->value;
            this->value = tmp.clone();
        }
        else
            (*this->value) = tmp;

        if ( this->dag_node != NULL )
        {
            this->dag_node->restoreAffected();
        }
    }
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::setValue(mixtureType *v, bool force)
{
    
    const RbVector<mixtureType> &vals = parameter_values->getValue();
    // we need to catch the value and increment the index
    for (index = 0; index < vals.size(); ++index)
    {
        if ( vals[index] == *v )
        {
            break;
        }
    }
    
    // delegate class
    TypedDistribution<mixtureType>::setValue( v, force );
}


template <class mixtureType>
void RevBayesCore::MixtureDistribution<mixtureType>::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    // only do this when the toucher was our parameters
    if ( toucher == parameter_values )
    {
        const mixtureType &tmp = parameter_values->getValue()[index];
        if constexpr (std::is_base_of_v<Cloneable, mixtureType>)
        {
            delete this->value;
            this->value = tmp.clone();
        }
        else
            (*this->value) = tmp;

        if ( this->dag_node != NULL )
        {
            this->dag_node->touchAffected();
        }
    }
    
}

#endif
