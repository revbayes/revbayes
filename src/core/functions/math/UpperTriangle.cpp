#include "UpperTriangle.h"
#include "RbException.h"

using namespace RevBayesCore;

UpperTriangle::UpperTriangle(const TypedDagNode< MatrixReal > *m) : TypedFunction< RbVector<double> >( new RbVector<double>(m->getValue().getDim(), 1.0) ),
    matrix( m )
{

    // check dimensions here
    addParameter( matrix );
    update();
    
}


UpperTriangle::~UpperTriangle( void ) {

}


UpperTriangle* UpperTriangle::clone( void ) const {

    return new UpperTriangle( *this );
}


void UpperTriangle::update(void) {

    const MatrixReal R = matrix->getValue();
    size_t nrows = R.getDim();

    value->clear();
    value->resize( nrows * (nrows - 1) / 2 );
    
    size_t k = 0;
    for (size_t i = 0; i < nrows; ++i)
    {
        for (size_t j = i + 1; j < nrows; ++j)
        {
            (*value)[k++] = R[i][j];
        }
    }

}


void UpperTriangle::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    // check dimensions here
    if (oldP == matrix)
    {
        matrix = static_cast<const TypedDagNode< MatrixReal >* >( newP );
    }
}


