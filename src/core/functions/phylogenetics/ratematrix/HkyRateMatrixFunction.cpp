#include "HkyRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


HkyRateMatrixFunction::HkyRateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_HKY() ),
    base_frequencies( bf ),
    kappa( k )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    addParameter( kappa );
    
    update();
}


HkyRateMatrixFunction::~HkyRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



HkyRateMatrixFunction* HkyRateMatrixFunction::clone( void ) const
{

    return new HkyRateMatrixFunction( *this );
}


void HkyRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    double k = kappa->getValue();
    const std::vector<double>& f = base_frequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_HKY* >(value)->setStationaryFrequencies( f );
    static_cast< RateMatrix_HKY* >(value)->setKappa( k );
    
    value->update();
}



void HkyRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}


