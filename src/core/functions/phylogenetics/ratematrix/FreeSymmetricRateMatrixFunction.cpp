#include "FreeSymmetricRateMatrixFunction.h"

#include <cstddef>
#include <cmath>
#include <vector>

#include "RateMatrix_FreeSymmetric.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Default constructor.
 *
 * This function takes three inputs:
 * @param tr The vector of substitution rates, whose number of elements equals half the number of the off-diagonal entries of the rate matrix.
 * @param r Should the rates be rescaled so that the average rate equals 1?
 * @param method Matrix exponentiation method.
 */

FreeSymmetricRateMatrixFunction::FreeSymmetricRateMatrixFunction(const TypedDagNode< RbVector<double> > *tr, bool r, std::string method) : TypedFunction<RateGenerator>( NULL ),
    transition_rates( tr )
{
    double num_off_diag_states = 2*tr->getValue().size();
    double tmp_num_states = 0.5 + sqrt(0.25+num_off_diag_states);
    size_t num_states = size_t( round( tmp_num_states ) );
    value = new RateMatrix_FreeSymmetric( num_states, r, method );
    
    // add the rate and frequency parameters as parents
    addParameter( transition_rates );
    
    update();
}


FreeSymmetricRateMatrixFunction::~FreeSymmetricRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



FreeSymmetricRateMatrixFunction* FreeSymmetricRateMatrixFunction::clone( void ) const
{
    return new FreeSymmetricRateMatrixFunction( *this );
}


void FreeSymmetricRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    const std::vector<double>& trr = transition_rates->getValue();
    
    // set the base frequencies
    static_cast< RateMatrix_FreeSymmetric* >(value)->setTransitionRates(trr);
    
    value->update();
}



void FreeSymmetricRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == transition_rates)
    {
        transition_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

