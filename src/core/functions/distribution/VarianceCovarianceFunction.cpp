#include "VarianceCovarianceFunction.h"

#include <cstddef>
#include <vector>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

VarianceCovarianceFunction::VarianceCovarianceFunction(const TypedDagNode< RbVector<double> > *sd, const TypedDagNode< RbVector<double> > *cc) : TypedFunction<MatrixReal>( new MatrixReal(sd->getValue().size()) ),
    standardDeviations( sd ),
    correlationCoefficients( cc )
{

    // check dimensions here
    // add the lambda parameter as a parent
    addParameter( standardDeviations );
    addParameter( correlationCoefficients );
    
    update();
}


VarianceCovarianceFunction::~VarianceCovarianceFunction( void ) {

}



VarianceCovarianceFunction* VarianceCovarianceFunction::clone( void ) const {

    return new VarianceCovarianceFunction( *this );
}


void VarianceCovarianceFunction::update(void)
{

    // get the information from the arguments for reading the file
    const std::vector<double>& s = standardDeviations->getValue();
    const std::vector<double>& p = correlationCoefficients->getValue();
    
    size_t k = 0;
    for (size_t i = 0; i < s.size(); i++)
    {
        (*value)[i][i] = s[i] * s[i];
        for (size_t j = i + 1; j < s.size(); j++)
        {
            (*value)[i][j] = s[i] * s[j] * p[k];
            (*value)[j][i] = s[i] * s[j] * p[k];
            k++;
        }
    }
}


void VarianceCovarianceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    // check dimensions here
    if (oldP == standardDeviations)
    {
        standardDeviations = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == correlationCoefficients)
    {
        correlationCoefficients = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}


