#include <cstddef>
#include <vector>

#include "MatrixReal.h"
#include "StochasticMatrixFunction.h"
#include "TypedDagNode.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** Standard constructor from appropriately typed DAG node */
StochasticMatrixFunction::StochasticMatrixFunction(const TypedDagNode<RbVector<Simplex > >* &args) :
TypedFunction< MatrixReal >( new MatrixReal() ),
matrixParams( args )
{

    // add the vector parameter as a parent
    this->addParameter( args );
    
    // initialize value
    const RbVector<Simplex >& v = matrixParams->getValue();
    *value = MatrixReal( v.size(), v[0].size(), 0.0 );
    
    // update the value
    update();
}


StochasticMatrixFunction::~StochasticMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



StochasticMatrixFunction* StochasticMatrixFunction::clone( void ) const
{
    return new StochasticMatrixFunction( *this );
}


void StochasticMatrixFunction::update( void )
{
	const RbVector<Simplex >& v = matrixParams->getValue();
	for (size_t i = 0; i < v.size(); i++) {
		for (size_t j = 0; j < v[i].size(); j++) {
			(*value)[i][j] = v[i][j];
		}
	}
}



void StochasticMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
	if (oldP == matrixParams)
	{
		matrixParams = static_cast<const TypedDagNode<RbVector<Simplex > >* >( newP );
	}
}
