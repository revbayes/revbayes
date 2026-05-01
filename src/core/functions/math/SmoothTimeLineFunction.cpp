#include "SmoothTimeLineFunction.h"

#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

SmoothTimeLineFunction::SmoothTimeLineFunction(const TypedDagNode< double > *max_t, const TypedDagNode< RbVector<double> > *t, const TypedDagNode< RbVector<double> > *v) : TypedFunction< RbVector<double> >( new RbVector<double>(v->getValue().size(),0.0) ),
    max_time( max_t ),
    times( t ),
    values( v )
{
    // add the lambda parameter as a parent
    addParameter( max_time );
    addParameter( times );
    addParameter( values );

    update();
}

SmoothTimeLineFunction::~SmoothTimeLineFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SmoothTimeLineFunction* SmoothTimeLineFunction::clone( void ) const
{
    return new SmoothTimeLineFunction( *this );
}


void SmoothTimeLineFunction::update( void )
{

    RbVector<double> &smooth_values = *value;
 
    const RbVector<double>& raw_values  = values->getValue();
    const RbVector<double>& time_values = times->getValue();
    double max_time_value = max_time->getValue();
    
    size_t num_times = time_values.size();
    
    // make sure we have the correct size
    if ( smooth_values.size() != num_times+1 )
    {
        smooth_values.resize( num_times + 1 );
    }
    
    smooth_values[0] = raw_values[0];

    // Assemble the field
    for (size_t i=1; i<num_times; ++i)
    {
        
        if ( max_time_value > time_values[i-1] )
        {
            smooth_values[i] = raw_values[i];
        }
        else
        {
            smooth_values[i] = smooth_values[i-1];
        }
        
    }

}

void SmoothTimeLineFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == times)
    {
        times = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == values)
    {
        values = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == max_time)
    {
        max_time = static_cast<const TypedDagNode< double >* >( newP );
    }

}
