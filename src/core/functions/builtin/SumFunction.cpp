#include "SumFunction.h"

#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "DagNode.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;
/**
 * SumFunction of a RbVector Constructor.
 * @param v the vector of values of type double
 */
SumFunction::SumFunction(const TypedDagNode<RbVector<double> > *v) : TypedFunction<double>( new double(0.0) ),
    matrix(false),
    vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


/** SumFunction of a MatrixReal Constructor
 * @param v the matrix of values
 */
SumFunction::SumFunction(const TypedDagNode< MatrixReal > *v) : TypedFunction<double>( new double(0.0) ),
    matrix(true),
    vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );

    update();
}


SumFunction::~SumFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SumFunction* SumFunction::clone( void ) const
{
    return new SumFunction( *this );
}


void SumFunction::update( void )
{
    
    double m = 0;
    if (matrix == true)
    {
        const MatrixReal &v = dynamic_cast<const TypedDagNode< MatrixReal >* >(vals)->getValue();
        
        for (size_t row = 0; row < v.size(); row++)
        {
            for (size_t col = 0; col < v[row].size(); col++)
            {
                m += v[row][col];
            }
        }
    }
    else
    {
        const RbVector<double> &v = dynamic_cast<const TypedDagNode< RbVector<double> >* >(vals)->getValue();
        
        for ( RbConstIterator<double> it = v.begin(); it != v.end(); ++it)
        {
            m += *it;
        }
    }
    
    *value = m ;
    
}



void SumFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = newP;
    }
    
}

