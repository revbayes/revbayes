#include "F81RateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


F81RateMatrixFunction::F81RateMatrixFunction(const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_F81(bf->getValue().size()) ),
    base_frequencies( bf )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    
    update();
}


F81RateMatrixFunction::~F81RateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



F81RateMatrixFunction* F81RateMatrixFunction::clone( void ) const
{
    return new F81RateMatrixFunction( *this );
}


void F81RateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    const std::vector<double>& f = base_frequencies->getValue();
    
    // set the base frequencies
    static_cast< RateMatrix_F81* >(value)->setStationaryFrequencies( f );
    
    value->update();
}



void F81RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
}


