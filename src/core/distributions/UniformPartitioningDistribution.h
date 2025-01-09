#ifndef UniformPartitioningDistribution_H
#define UniformPartitioningDistribution_H

#include "MemberObject.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
    /**
     * The uniform partitioning distribution samples uniformly distributed random partitions
     * of a set.
     *
     * The input is a vector of non-zero values that represents the most finely fragmented
     * partition of the set. There are B_n possible partitions, where n is the number of 
     * elements and B_n is the Bell number for n. Each of the B_n partitions have equal 
     * probability.
     *
     * Additionally, zero values can be added, so we now include all partitions of n elements
     * where each subset can either have one of the n non-zero values or a zero value. This means 
     * there are now B_{n+1} possible partitions, each of which have equal probability.
     *
     * Will Freyman 6/29/17
     *
     */
    template <class valueType>
    class UniformPartitioningDistribution : public TypedDistribution< RbVector<valueType> >, public MemberObject<long> {
        
    public:
        // constructor(s)
        UniformPartitioningDistribution(const TypedDagNode< RbVector<valueType> > *v, bool include_zero);
        UniformPartitioningDistribution(const UniformPartitioningDistribution &d);

        // public member functions
        UniformPartitioningDistribution*                    clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, long &rv) const;     //!< Map the member methods to internal function calls
        int                                                 getCurrentIndex(void) const;
        RbVector<valueType>*                                getPartition(void);
        std::vector<int>                                    getValueAssignments(void);

        void                                                redrawValue(void);
        void                                                setCurrentIndex(int i);
        void                                                setValue(valueType *v, bool f=false);
        void                                                setValueAssignments(std::vector<int> v);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        // helper methods
        RbVector<valueType>*                                simulate();
        RbVector<valueType>*                                assignValues();
        
        // private members
        bool                                                include_zero;
        int                                                 index;
        int                                                 num_partitions;
        const TypedDagNode< RbVector<valueType> >*          parameter_values;
        std::vector<int>                                    partition_label;
        std::vector<int>                                    value_assignments;
    };
    
}

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"

#include <cmath>
#include <vector>


template <class valueType>
RevBayesCore::UniformPartitioningDistribution<valueType>::UniformPartitioningDistribution(const TypedDagNode< RbVector<valueType> > *v, bool include_zero) : TypedDistribution< RbVector<valueType> >( new RbVector<valueType>() ),
    include_zero( include_zero ),
    index( 0 ),
    parameter_values( v )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( parameter_values );
    
    // calculate the number of partitions
    int num_elements = (int)parameter_values->getValue().size();
    if ( include_zero == false )
    {
        num_partitions = RbMath::bell( num_elements ); 
    }
    else
    {
        num_partitions = RbMath::bell( num_elements + 1 );
    }

    value_assignments = std::vector<int>(num_elements, 0);
    partition_label = std::vector<int>(num_elements, 0);
    
    delete this->value;
    this->value = simulate();
}


template <class valueType>
RevBayesCore::UniformPartitioningDistribution<valueType>::UniformPartitioningDistribution( const UniformPartitioningDistribution &d ) : TypedDistribution< RbVector<valueType> >( d ),
    include_zero( d.include_zero ),
    index( d.index ),
    parameter_values( d.parameter_values ),
    partition_label( d.partition_label ),
    value_assignments( d.value_assignments )
{

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( parameter_values );
    
    // calculate the number of partitions
    int num_elements = (int)parameter_values->getValue().size();
    if ( include_zero == false )
    {
        num_partitions = RbMath::bell( num_elements ); 
    }
    else
    {
        num_partitions = RbMath::bell( num_elements + 1 );
    }
    
    delete this->value;
    this->value = simulate();

}

    
template <class valueType>
RevBayesCore::RbVector<valueType>* RevBayesCore::UniformPartitioningDistribution<valueType>::assignValues( void )
{
    
    int num_elements = (int)parameter_values->getValue().size();
    RbVector<valueType>* rv = new RbVector<valueType>();
    RbVector<valueType>& temp = *rv;
    
    if ( include_zero == true )
    {
        for (int i = 1; i < (num_elements + 1); i++)
        {
            if (partition_label[i] == 0)
                temp.push_back(0.0);
            else
            {
                valueType val = parameter_values->getValue()[ value_assignments[ partition_label[i] - 1 ] ];
                temp.push_back( val );
            }
        }
    }
    else
    {
        for (int i = 0; i < num_elements; i++)
        {
            valueType val = parameter_values->getValue()[ value_assignments[ partition_label[i] - 1 ] ];
            temp.push_back( val );
        }
    }

    return rv;

}


template <class valueType>
RevBayesCore::UniformPartitioningDistribution<valueType>* RevBayesCore::UniformPartitioningDistribution<valueType>::clone( void ) const
{

    return new UniformPartitioningDistribution<valueType>( *this );
}



template <class valueType>
double RevBayesCore::UniformPartitioningDistribution<valueType>::computeLnProbability( void )
{
    
    return -log( num_partitions );

}


template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, long &rv) const
{
    
    if ( n == "getAllocationIndex" )
    {
        rv = long(index) + 1;
    }
    else
    {
        throw RbException() << "A uniform partitions distribution does not have a member method called '" << n << "'.";
    }
    
}


template <class valueType>
int RevBayesCore::UniformPartitioningDistribution<valueType>::getCurrentIndex( void ) const
{

    return index;
}


template <class valueType>
RevBayesCore::RbVector<valueType>* RevBayesCore::UniformPartitioningDistribution<valueType>::getPartition( void )
{
    int num_elements = (int)parameter_values->getValue().size();

    if ( include_zero == true )
    {
        partition_label = std::vector<int>(num_elements + 1, 0);
    }
    else
    {
        partition_label = std::vector<int>(num_elements, 1);
    }

    // iterate through an ordered list of all partitions, stopping on the current index
    // TODO list of partitions should be generated just the first time and then cached....
    int current_partition = 0;
    while (current_partition < int(index))
    //while (current_partition < num_partitions)
    {
        int start_index;
        if ( include_zero == true )
        {
            start_index = num_elements;
        }
        else
        {
            start_index = num_elements - 1;
        }

        while (start_index >= 0)
        {
            int max_value = 0;
            for (int i = 0; i < start_index; i++)
            {
                if (partition_label[i] > max_value)
                    max_value = partition_label[i];
            }

            if (partition_label[start_index] > max_value or max_value > num_elements or partition_label[start_index] >= num_elements)
            {
                if ( include_zero == true )
                {
                    partition_label[start_index] = 0;
                }
                else
                {
                    partition_label[start_index] = 1;
                }
                start_index--;
            }
            else
            {
                partition_label[start_index] += 1;
                break;
            }
        }

        current_partition++;
    }
   
    // now we must randomly assign values to the partition's non-zero subsets
    
    // first find the max element label
    int max_value = 0;
    for (int i = 0; i < partition_label.size(); i++)
    {
        if (partition_label[i] > max_value)
            max_value = partition_label[i];
    }
    
    // for each subset label (1,2,3, etc) assign one of the non-zero elements
    // value_assignments[ label - 1 ] = element_num
    double weight = 1 / double(num_elements);
    RandomNumberGenerator *rng = GLOBAL_RNG;
    value_assignments = std::vector<int>();
    for (int i = 0; i < max_value; i++)
    {
        while (true)
        {
            double u = rng->uniform01();
            int j = 0;
            while ( u > weight )
            {
                u -= weight;
                ++j;
            }
            bool already_used = false;
            for (int k = 0; k < value_assignments.size(); k++)
            {
                if (value_assignments[k] == j)
                {
                    already_used = true;
                    break;
                }
            }
            if (already_used == false)
            {
                value_assignments.push_back(j);
                break;
            }
        }
    }

    return assignValues();

}


template <class valueType>
std::vector<int> RevBayesCore::UniformPartitioningDistribution<valueType>::getValueAssignments( void )
{
    return value_assignments;
}


template <class valueType>
RevBayesCore::RbVector<valueType>* RevBayesCore::UniformPartitioningDistribution<valueType>::simulate()
{
    // draw a random partition index 
    double weight = 1 / double(num_partitions);

    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    
    index = 0;
    while ( u > weight )
    {
        u -= weight;
        ++index;
    }
    
    // now get that partition
    return getPartition();
}


template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::setCurrentIndex( int i )
{
    index = i;
    delete this->value;
    this->value = getPartition();
}


template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::setValueAssignments(std::vector<int> v)
{
    value_assignments = v;
    delete this->value;
    this->value = assignValues();
}


template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::redrawValue( void )
{
    delete this->value;
    this->value = simulate();
    //Assign<valueType, IsDerivedFrom<valueType, Assignable>::Is >::doAssign( (*this->value), simulate() );
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == parameter_values)
    {
        parameter_values = static_cast<const TypedDagNode< RbVector<valueType> >* >( newP );
    }
}


template <class valueType>
void RevBayesCore::UniformPartitioningDistribution<valueType>::setValue(valueType *v, bool force)
{
    
    const RbVector<valueType> &vals = parameter_values->getValue();
    // we need to catch the value and increment the index
    for (index = 0; index < vals.size(); ++index)
    {
        if ( vals[index] == *v )
        {
            break;
        }
    }
    
    // delegate class
    TypedDistribution<valueType>::setValue( v, force );
}



#endif
