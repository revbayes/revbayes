#include "AvgDistanceMatrixFunction.h"

#include <vector>

#include "RbVector.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

AvgDistanceMatrixFunction::AvgDistanceMatrixFunction(const TypedDagNode< RbVector<DistanceMatrix> >* matvect) : TypedFunction< AverageDistanceMatrix >( new AverageDistanceMatrix() ),
matrixVector( matvect )
{
    weightVector = NULL;
    addParameter( matrixVector );
    
    update();
}


AvgDistanceMatrixFunction::AvgDistanceMatrixFunction(const TypedDagNode< RbVector<DistanceMatrix> >* matvect, const TypedDagNode< RbVector<double> >* weights) : TypedFunction< AverageDistanceMatrix >( new AverageDistanceMatrix() ),
matrixVector( matvect ),
weightVector( weights )
{
    addParameter( matrixVector );
    addParameter( weightVector );
    
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
    
    if (weightVector == NULL)
    {
        *value = TreeUtilities::getAverageDistanceMatrix( matrixVector->getValue(), NULL );
    }
    else
    {
        *value = TreeUtilities::getAverageDistanceMatrix( matrixVector->getValue(), &weightVector->getValue() );
    }
}


void AvgDistanceMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == matrixVector)
    {
        matrixVector = static_cast<const TypedDagNode< RbVector<DistanceMatrix> >* >( newP );
    }
    else if (oldP == weightVector)
    {
        weightVector = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
