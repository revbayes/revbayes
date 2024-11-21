#ifndef AnalyticalMixtureDistribution_H
#define AnalyticalMixtureDistribution_H

#include "MemberObject.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements a generic mixture distribution between several possible values.
     *
     * This mixture can be considered as a multinomial distribution. We specify a Analytical of probabilities
     * and a Analytical of values. Then, a value drawn from this distribution takes each value corresponding to
     * its probability.
     * The values are already of the correct mixture type. You may want to apply a mixture allocation move
     * to change between the current value. The values themselves change automatically when the input parameters change.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    template <class mixtureType>
    class AnalyticalMixtureDistribution : public TypedDistribution< mixtureType >, public MemberObject< Simplex > {
        
    public:
        // constructor(s)
        AnalyticalMixtureDistribution(std::vector<TypedDistribution< mixtureType > *> base_dist, const TypedDagNode< Simplex > *p);
        AnalyticalMixtureDistribution(const AnalyticalMixtureDistribution<mixtureType> &d);
        
        // public member functions
        AnalyticalMixtureDistribution*                          clone(void) const;                                                                      //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Simplex &rv) const;     //!< Map the member methods to internal function calls
        void                                                    redrawValue(void);
        //        void                                                    setValue(RbAnalytical<mixtureType> *v, bool f=false);
        
        // special handling of state changes
        void                                                    getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                          //!< get affected nodes
        void                                                    keepSpecialization(const DagNode* affecter);
        void                                                    restoreSpecialization(const DagNode *restorer);
        void                                                    touchSpecialization(const DagNode *toucher, bool touchAll);
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        // helper methods
        
        // private members
        std::vector<TypedDistribution< mixtureType >* >         base_distributions;
        const TypedDagNode< Simplex >*                          probabilities;
        
        bool                                                    dirty;
        std::vector<double>                                     ln_probabilities;
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>

template <class mixtureType>
RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::AnalyticalMixtureDistribution(std::vector<TypedDistribution< mixtureType > *> base_dists, const TypedDagNode< Simplex > *p) : TypedDistribution< mixtureType >( new mixtureType() ),
    base_distributions( base_dists ),
    probabilities( p )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( probabilities );
    
    // add the parameters of the distribution
    for ( size_t i=0; i<base_distributions.size(); ++i)
    {
        const std::vector<const DagNode*>& pars = base_distributions[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }
    
    ln_probabilities = std::vector<double>(base_distributions.size(), 0.0);
    dirty = true;
    
    this->redrawValue();
}


template <class mixtureType>
RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::AnalyticalMixtureDistribution( const AnalyticalMixtureDistribution<mixtureType> &d ) : TypedDistribution< mixtureType >(d),
    base_distributions(),
    probabilities( d.probabilities ),
    dirty( d.dirty ),
    ln_probabilities( d.ln_probabilities )
{
    
    // add the parameters of the distribution
    for ( size_t i=0; i<d.base_distributions.size(); ++i )
    {
        // first we need to clone the base distribution
        base_distributions.push_back( d.base_distributions[i]->clone() );
        
        const std::vector<const DagNode*>& prior_pars = base_distributions[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = prior_pars.begin(); it != prior_pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }
    
}



template <class mixtureType>
RevBayesCore::AnalyticalMixtureDistribution<mixtureType>* RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::clone( void ) const
{
    
    return new AnalyticalMixtureDistribution<mixtureType>( *this );
}



template <class mixtureType>
double RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::computeLnProbability( void )
{
    dirty = false;
    
    const Simplex &probs = probabilities->getValue();
    for ( size_t i=0; i<base_distributions.size(); ++i )
    {
        // get the i-th base distribution
        TypedDistribution<mixtureType> *this_base_dist = base_distributions[i];
        
        // multiply with the prior probability for this base distribution (i.e., mixture type)
        ln_probabilities[i] = log( probs[i] );
        
        // now compute the probability for each value under this base distribution
        const mixtureType &tmp_val = (*this->value);
        this_base_dist->setValue( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( tmp_val ) );
        ln_probabilities[i] = this_base_dist->computeLnProbability();

    }
    
    double max_prob = ln_probabilities[0];
    for (size_t i=1; i<base_distributions.size(); ++i)
    {
        if ( ln_probabilities[i] > max_prob )
        {
            max_prob = ln_probabilities[i];
        }
    }
    
    // safe summation of exp-transformed values
    double sum_probs = 0.0;
    for (size_t i=0; i<base_distributions.size(); ++i)
    {
        sum_probs += exp( ln_probabilities[i] - max_prob );
    }
    double ln_sum_probs = log( sum_probs ) + max_prob;
    
    
    return ln_sum_probs;
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, Simplex &rv) const
{
    
    if ( n == "getMixtureProbabilities" )
    {
        rv = Simplex( base_distributions.size() );
        
        if (dirty == true)
        {
            const_cast< AnalyticalMixtureDistribution<mixtureType>* >(this)->computeLnProbability();
        }
        
        // first, find the maximum probability
        double max_prob = ln_probabilities[0];
        for (size_t i=1; i<base_distributions.size(); ++i)
        {
            if ( ln_probabilities[i] > max_prob )
            {
                max_prob = ln_probabilities[i];
            }
        }
        
        // next, compute the total probability needed for normalization
        double sum_probs = 0.0;
        for (size_t i=0; i<base_distributions.size(); ++i)
        {
            sum_probs += exp( ln_probabilities[i] - max_prob );
        }
        double ln_sum_probs = log( sum_probs ) + max_prob;
        
        // finally, compute each probability
        for (size_t i=0; i<base_distributions.size(); ++i)
        {
            rv[i] = exp( ln_probabilities[i] - ln_sum_probs );
        }
        
    }
    else
    {
        throw RbException() << "A Analytical-mixture distribution does not have a member method called '" << n << "'.";
    }
    
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter)
{
    
    
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::keepSpecialization( const DagNode* affecter )
{
    // only do this when the toucher was our parameters
    //    if ( affecter == parameterValues && this->dag_node != NULL )
    //    {
    //        this->dag_node->keepAffected();
    //    }
    
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::redrawValue( void )
{
    
    const Simplex &probs = probabilities->getValue();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    size_t index = 0;
    while ( u > probs[index] )
    {
        u -= probs[index];
        ++index;
    }
    
    TypedDistribution<mixtureType> *selected_base_dist = base_distributions[index];
    selected_base_dist->redrawValue();
    if constexpr (std::is_base_of_v<Cloneable, mixtureType>)
    {
	delete this->value;
	this->value = selected_base_dist->getValue().clone();
    }
    else
	(*this->value) = selected_base_dist->getValue();
}


/** Swap a parameter of the distribution */
template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::swapParameterInternal( const DagNode *old_p, const DagNode *new_p )
{
    bool found = false;
    
    if (old_p == probabilities)
    {
        found = true;
        probabilities = static_cast<const TypedDagNode< Simplex >* >( new_p );
    }
    else
    {
        for (int j = 0; j < base_distributions.size(); ++j)
        {
            
            TypedDistribution<mixtureType> *this_prior = base_distributions[j];
            try {
                this_prior->swapParameter(old_p,new_p);
                // if the statement succeeded and didn't throw an error, then the distribution had this parameter
                found = true;
            }
            catch (RbException &e)
            {
                // do nothing because we actually do not know who had the parameter
            }
            
        }
        
    }
    
    if ( found == false )
    {
        throw RbException("Could not find the distribution parameter to be swapped: " + old_p->getName() + " to " + new_p->getName()) ;
    }
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::restoreSpecialization( const DagNode *restorer )
{
    
    // only do this when the toucher was our parameters
    dirty = true;
    
    
}


template <class mixtureType>
void RevBayesCore::AnalyticalMixtureDistribution<mixtureType>::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    // only do this when the toucher was our parameters
    dirty = true;
    
}

#endif
