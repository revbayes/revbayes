#include "DecomposedVarianceCovarianceFunction.h"

#include <cstddef>
#include <vector>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

DecomposedVarianceCovarianceFunction::DecomposedVarianceCovarianceFunction(const TypedDagNode< RbVector<double> > *sd, const TypedDagNode< MatrixReal > *cc) : TypedFunction<MatrixReal>( new MatrixReal(sd->getValue().size()) ),
    standardDeviations( sd ),
    correlationCoefficients( cc ) {

    // check dimensions here
    // add the lambda parameter as a parent
    addParameter( standardDeviations );
    addParameter( correlationCoefficients );
    
    update();
}


DecomposedVarianceCovarianceFunction::~DecomposedVarianceCovarianceFunction( void ) {

}


DecomposedVarianceCovarianceFunction* DecomposedVarianceCovarianceFunction::clone( void ) const {

    return new DecomposedVarianceCovarianceFunction( *this );
}


void DecomposedVarianceCovarianceFunction::update(void) {

    // get the information from the arguments for reading the file
    const std::vector<double>& sd = standardDeviations->getValue();
    const MatrixReal R = correlationCoefficients->getValue();
    size_t dim = sd.size();
    
    // Fill in the matrix with the standard deviations. The off diagonal elements
    // are zero while the diagonal elements contain the standard deviations.
    MatrixReal S = MatrixReal(dim);
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            {
            if (i == j)
                S[i][j] = sd[i];
            else
                S[i][j] = 0.0;
            }
        }
    
    // We want to do the following -- S * R * S -- where S is the matrix of
    // standard deviations and R is a matrix with the correlation coefficients
    // on the off-diagonal components and 1.0 as the diagonal components.
    // First, multiply S * R ...
    MatrixReal temp = MatrixReal(dim);
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            {
            double sum = 0.0;
            for (int k=0; k<dim; k++)
                sum += S[i][k] * R[k][j];
            temp[i][j] = sum;
            }
        }
    
    // ... then, multiply that result by (S*R) * S:
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            {
            double sum = 0.0;
            for (int k=0; k<dim; k++)
                sum += temp[i][k] * S[k][j];
            (*value)[i][j] = sum;
            }
        }
}


void DecomposedVarianceCovarianceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    // check dimensions here
    if (oldP == standardDeviations)
        {
        standardDeviations = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
        }
    else if (oldP == correlationCoefficients)
        {
        correlationCoefficients = static_cast<const TypedDagNode< MatrixReal >* >( newP );
        }
}


