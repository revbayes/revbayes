#include <cstddef>
#include <cmath>

#include "PartialToCorrelationFunction.h"
#include "Cloner.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PartialToCorrelation::PartialToCorrelation(const TypedDagNode< MatrixReal > *pc) : TypedFunction<MatrixReal>( new MatrixReal(pc->getValue().getNumberOfRows()) ),
    partial_correlations( pc )
{

    // check dimensions here
    // add the lambda parameter as a parent
    addParameter( partial_correlations );
    
    update();
}


PartialToCorrelation::~PartialToCorrelation( void ) {

}


PartialToCorrelation* PartialToCorrelation::clone( void ) const {

    return new PartialToCorrelation( *this );
}


void PartialToCorrelation::update(void) {

    // compute the correlation matrix
    const MatrixReal P = partial_correlations->getValue();
    size_t dim = P.getNumberOfRows();
    
    // P must start as an identity matrix
    for (size_t i = 0; i < dim; ++i)
    {
        (*value)[i][i] = 1.0;
    }
    
    for (int k = 0; k < dim - 1; ++k)
    {
        for (int i = k + 1; i < dim; ++i)
        {
            double p = P[k][i];
            for (int l = k - 1; l >= 0; --l)
            {
                p = p * pow( ( 1.0 - pow(P[l][i], 2) ) * ( 1.0 - pow(P[l][k], 2) ), 0.5) + P[l][i] * P[l][k];
            }
            (*value)[i][k] = p;
            (*value)[k][i] = p;
        }
    }
    
    return;
    
}


void PartialToCorrelation::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    // check dimensions here
    if (oldP == partial_correlations)
    {
        partial_correlations = static_cast<const TypedDagNode< MatrixReal >* >( newP );
    }

}


