#include "AvgDistanceMatrixFunction.h"

#include <vector>

#include "RbVector.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

AvgDistanceMatrixFunction::AvgDistanceMatrixFunction(const TypedDagNode< RbVector<DistanceMatrix> > *matvect) : TypedFunction< AverageDistanceMatrix >( new AverageDistanceMatrix() ),
    matrixVector( matvect )
{
    addParameter( matrixVector );
    
    update();
}


AvgDistanceMatrixFunction::~AvgDistanceMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


AvgDistanceMatrixFunction* AvgDistanceMatrixFunction::clone( void ) const
{
    
    return new AvgDistanceMatrixFunction( *this );
}


void AvgDistanceMatrixFunction::update( void )
{

    *value = TreeUtilities::getAverageDistanceMatrix( matrixVector->getValue() );
}


void AvgDistanceMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == matrixVector)
    {
        matrixVector = static_cast<const TypedDagNode< RbVector<DistanceMatrix> >* >( newP );
    }
    
}
