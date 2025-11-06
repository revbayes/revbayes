#include "NormalizeVectorFunction.h"

#include <cstddef>
#include <vector>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

NormalizeVectorFunction::NormalizeVectorFunction(const TypedDagNode< RbVector<double> > *v, const TypedDagNode< double > *s) : TypedFunction< RbVector<double> >( new RbVector<double>(v->getValue().size(), 1.0) ),
    vals( v ),
    sum( s )
{
    // add the parameters as parents
    this->addParameter( vals );
    this->addParameter( sum );
    
    update();
}


NormalizeVectorFunction* NormalizeVectorFunction::clone( void ) const
{   
    return new NormalizeVectorFunction( *this );
}


void NormalizeVectorFunction::update( void )
{
    
    double m = 0;
    const std::vector<double> &v = vals->getValue();
    for ( double val : v)
    {
        m += val;
    }
    
    m /= v.size() * sum->getValue();
    
    
    for ( size_t i = 0; i < v.size(); ++i ) 
    {
        (*value)[i] = v[i] / m;
    }
    
    
}



void NormalizeVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals ) 
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == sum )
    {
        sum = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}

