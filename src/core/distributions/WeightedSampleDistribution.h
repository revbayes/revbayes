#ifndef WeightedSampleDistribution_H
#define WeightedSampleDistribution_H

#include "MemberObject.h"
#include "Parallelizable.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements a generic weighted-sample distribution.
     *
     * This distribution represents a wrapper distribution for basically any other distribution.
     * This distribution should be used when the "observed" value is not known exactly but instead
     * is known by a sample from, e.g., its posterior distribution. Then, we compute the probability
     * of each value and compute the mean probability.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    template <class valueType>
    class WeightedSampleDistribution : public TypedDistribution< RbVector<valueType> > {
        
    public:
        // constructor(s)
        WeightedSampleDistribution(TypedDistribution<valueType> *g, std::vector<TypedDistribution<valueType>* > g_vector, size_t ns);
        WeightedSampleDistribution(const WeightedSampleDistribution &d);
        virtual                                            ~WeightedSampleDistribution(void);
        
        WeightedSampleDistribution&                         operator=(const WeightedSampleDistribution &d);
        
        // public member functions
        WeightedSampleDistribution*                         clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        RevLanguage::RevPtr<RevLanguage::RevVariable>       executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found);
        void                                                redrawValue(void);
        void                                                setValue(RbVector<valueType> *v, bool f=false);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        // helper methods
        RbVector<valueType>*                                simulate();
        
        // private members
        TypedDistribution<valueType>*                       base_distribution;
//        std::vector< TypedDistribution<valueType>* >        base_distribution_vector;
        std::vector< TypedDistribution<valueType>* >        base_distribution_instances;
        
        bool                                                use_base_distribution_iid;
        
        size_t                                              num_samples;
        size_t                                              sample_block_start;
        size_t                                              sample_block_end;
        size_t                                              sample_block_size;
        
        
#ifdef RB_MPI
        std::vector<size_t>                                 pid_per_sample;
#endif
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <cmath>


#ifdef RB_MPI
#include <mpi.h>
#endif

template <class valueType>
RevBayesCore::WeightedSampleDistribution<valueType>::WeightedSampleDistribution(TypedDistribution<valueType> *g, std::vector<TypedDistribution<valueType>* > g_vector, size_t ns) : TypedDistribution< RbVector<valueType> >( new RbVector<valueType>() ),
base_distribution( g ),
base_distribution_instances( g_vector ),
use_base_distribution_iid(true),
num_samples( ns ),
sample_block_start( 0 ),
sample_block_end( num_samples ),
sample_block_size( num_samples )
{
    if (g_vector.size() > 0)
    {
        use_base_distribution_iid = false;
        base_distribution = NULL;
    }
    
    if (use_base_distribution_iid)
    {
        // add the parameters of the base distribution
        const std::vector<const DagNode*>& pars = base_distribution->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            const DagNode *the_node = *it;
            this->addParameter( the_node );
        }
    }
    else
    {
        // add the parameters for each base distribution
        for (size_t i = 0; i < num_samples; i++) {
            const std::vector<const DagNode*>& pars = base_distribution_instances[i]->getParameters();
            for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
            {
                const DagNode *the_node = *it;
                this->addParameter( the_node );
            }
        }
    }
    
    delete this->value;
    
    RbVector<valueType>* new_value = simulate();
    this->setValue( new_value );
    
}

template <class valueType>
RevBayesCore::WeightedSampleDistribution<valueType>::WeightedSampleDistribution( const WeightedSampleDistribution &d ) : TypedDistribution< RbVector<valueType> >( new RbVector<valueType>( d.getValue() ) ),
base_distribution( d.base_distribution->clone() ),
base_distribution_instances( ),
use_base_distribution_iid(d.use_base_distribution_iid),
num_samples( d.num_samples ),
sample_block_start( d.sample_block_start ),
sample_block_end( d.sample_block_end ),
sample_block_size( d.sample_block_size )
{
    
#ifdef RB_MPI
    pid_per_sample = d.pid_per_sample;
#endif
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        const DagNode *the_node = *it;
        this->addParameter( the_node );
    }
    
    base_distribution_instances = std::vector< TypedDistribution<valueType>* >( num_samples, NULL );
    for (size_t i = sample_block_start; i < sample_block_end; ++i)
    {
        base_distribution_instances[i] = d.base_distribution_instances[i]->clone();
    }
    
}


template <class valueType>
RevBayesCore::WeightedSampleDistribution<valueType>::~WeightedSampleDistribution( void )
{
    if (base_distribution != NULL)
    {
        delete base_distribution;
    }
    for (size_t i = 0; i < num_samples; ++i)
    {
        delete base_distribution_instances[i];
    }
    
}


template <class valueType>
RevBayesCore::WeightedSampleDistribution<valueType>& RevBayesCore::WeightedSampleDistribution<valueType>::operator=(const WeightedSampleDistribution &d)
{
    
    if ( this != &d )
    {
        // call base class
        TypedDistribution< RbVector<valueType> >::operator=( d );
        
        delete base_distribution;
        for (size_t i = 0; i < num_samples; ++i)
        {
            delete base_distribution_instances[i];
        }
        base_distribution_instances.clear();
        use_base_distribution_iid = d.use_base_distribution_iid;
        
        num_samples         = d.num_samples;
        sample_block_start  = d.sample_block_start;
        sample_block_end    = d.sample_block_end;
        sample_block_size   = d.sample_block_size;
        
#ifdef RB_MPI
        pid_per_sample = d.pid_per_sample;
#endif

        if (use_base_distribution_iid) {
            base_distribution   = d.base_distribution->clone();
        }
        
        // populate the base distribution instances
        base_distribution_instances = std::vector< TypedDistribution<valueType>* >( num_samples, NULL );
        for (size_t i = sample_block_start; i < sample_block_end; ++i)
        {
            base_distribution_instances[i] = d.base_distribution_instances[i]->clone();
        }
        
        if (use_base_distribution_iid)
        {
            // add the parameters of the base distribution
            const std::vector<const DagNode*>& pars = base_distribution->getParameters();
            for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
            {
                this->addParameter( *it );
            }
        }
        else
        {
            // add the parameters for each base distribution
            for (size_t i = 0; i < num_samples; i++) {
                const std::vector<const DagNode*>& pars = base_distribution_instances[i]->getParameters();
                for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
                {
                    this->addParameter( *it );
                }
            }
        }        
    }
    
    return *this;
}


template <class valueType>
RevBayesCore::WeightedSampleDistribution<valueType>* RevBayesCore::WeightedSampleDistribution<valueType>::clone( void ) const
{
    
    return new WeightedSampleDistribution<valueType>( *this );
}



template <class valueType>
double RevBayesCore::WeightedSampleDistribution<valueType>::computeLnProbability( void )
{
    
    double ln_prob = 0;
    double prob    = 0;
    
    std::vector<double> ln_probs = std::vector<double>(num_samples, 0.0);
    std::vector<double> probs    = std::vector<double>(num_samples, 0.0);
    
    // add the ln-probs for each sample
    for (size_t i = 0; i < num_samples; ++i)
    {
        if ( i >= sample_block_start && i < sample_block_end )
        {
            ln_probs[i] = base_distribution_instances[i]->computeLnProbability();
        }
        
    }
    
    
#ifdef RB_MPI
    for (size_t i = 0; i < num_samples; ++i)
    {
        
        if ( this->pid == pid_per_sample[i] )
        {
            
            // send the likelihood from the helpers to the master
            if ( this->process_active == false )
            {
                // send from the workers the log-likelihood to the master
                MPI_Send(&ln_probs[i], 1, MPI_DOUBLE, this->active_PID, 0, MPI_COMM_WORLD);
            }
            
        }
        // receive the likelihoods from the helpers
        else if ( this->process_active == true )
        {
            MPI_Status status;
            MPI_Recv(&ln_probs[i], 1, MPI_DOUBLE, pid_per_sample[i], 0, MPI_COMM_WORLD, &status);
        }
        
    }
#endif
    
    
    double max = 0;
    // add the ln-probs for each sample
    for (size_t i = 0; i < num_samples; ++i)
    {
        
#ifdef RB_MPI
        if ( this->process_active == true )
        {
#endif
            if ( i == 0 || max < ln_probs[i] )
            {
                max = ln_probs[i];
            }
            
#ifdef RB_MPI
        }
#endif
        
    }
    
    
#ifdef RB_MPI
    if ( this->process_active == true )
    {
#endif
        
        // now normalize
        for (size_t i = 0; i < num_samples; ++i)
        {
            probs[i] = exp( ln_probs[i] - max);
            prob += probs[i];
        }
        
        ln_prob = std::log( prob ) + max - std::log( num_samples );
        
#ifdef RB_MPI
        
        for (size_t i=this->active_PID+1; i<this->active_PID+this->num_processes; ++i)
        {
            MPI_Send(&ln_prob, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD);
        }
        
    }
    else
    {
        
        MPI_Status status;
        MPI_Recv(&ln_prob, 1, MPI_DOUBLE, this->active_PID, 0, MPI_COMM_WORLD, &status);
        
    }
#endif
    
    return ln_prob;
}


template <class valueType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevBayesCore::WeightedSampleDistribution<valueType>::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
{
    
    bool org_found = found;
    for (size_t i = 0; i < num_samples; ++i)
    {
        bool f = org_found;
        if ( base_distribution_instances[i] != NULL )
        {
            base_distribution_instances[i]->executeProcedure(name, args, f);
        }
        found |= f;
    }
    return NULL; 
}


template <class valueType>
RevBayesCore::RbVector<valueType>* RevBayesCore::WeightedSampleDistribution<valueType>::simulate()
{
    
    RbVector<valueType> *values = new RbVector<valueType>( num_samples );
    for (size_t i = 0; i < num_samples; ++i)
    {
        base_distribution_instances[i]->redrawValue();
        (*values)[i] = base_distribution_instances[i]->getValue();
    }
    
    return values;
}


template <class valueType>
void RevBayesCore::WeightedSampleDistribution<valueType>::redrawValue( void )
{
    
    delete this->value;
    RbVector<valueType> *new_value = simulate();
    this->setValue( new_value );
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::WeightedSampleDistribution<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (use_base_distribution_iid)
    {
        base_distribution->swapParameter(oldP,newP);
    }
    
    for (size_t i = 0; i < num_samples; ++i)
    {
        if ( base_distribution_instances[i] != NULL )
        {
            base_distribution_instances[i]->swapParameter(oldP,newP);
        }
        
    }
    
}


template <class valueType>
void RevBayesCore::WeightedSampleDistribution<valueType>::setValue(RbVector<valueType> *v, bool force)
{
    
    // free the old distributions
    
    num_samples = v->size();
    
    // compute which block of the data this process needs to compute
    sample_block_start = 0;
    sample_block_end   = num_samples;
#ifdef RB_MPI
    sample_block_start = size_t(floor( (double(this->pid-this->active_PID)   / this->num_processes ) * num_samples) );
    sample_block_end   = size_t(floor( (double(this->pid+1-this->active_PID) / this->num_processes ) * num_samples) );
#endif
    sample_block_size  = sample_block_end - sample_block_start;
    
    if (use_base_distribution_iid)
    {
        for (size_t i = 0; i < num_samples; ++i)
        {
            delete base_distribution_instances[i];
        }

        base_distribution_instances = std::vector< TypedDistribution<valueType>* >( num_samples, NULL );
        for (size_t i = sample_block_start; i < sample_block_end; ++i)
        {
            TypedDistribution<valueType> *base_distribution_clone = base_distribution->clone();
            base_distribution_instances[i] = base_distribution_clone;
            base_distribution_clone->setValue( (*v)[i].clone() );
        }
    }
    else
    {
        for (size_t i = sample_block_start; i < sample_block_end; ++i)
        {
//            TypedDistribution<valueType> *base_distribution_clone = base_distribution->clone();
//            base_distribution_instances[i] = base_distribution_clone;
            base_distribution_instances[i]->setValue( (*v)[i].clone() );
        }
    }
    
#ifdef RB_MPI
    // now we need to populate how is
    pid_per_sample = std::vector<size_t>( num_samples, 0 );
    for (size_t i = 0; i < num_samples; ++i)
    {
        pid_per_sample[i] = size_t(floor( double(i) / num_samples * this->num_processes ) ) + this->active_PID;
    }
#endif
    
    // delegate class
    TypedDistribution< RbVector<valueType> >::setValue( v, force );
    
}

#endif
