#include "K80RateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

K80RateMatrixFunction::K80RateMatrixFunction(const TypedDagNode<double> *k) : TypedFunction<RateGenerator>( new RateMatrix_Kimura80() ),
    kappa( k )
{
    // add the lambda parameter as a parent
    addParameter( kappa );
    
    update();
}


K80RateMatrixFunction::~K80RateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



K80RateMatrixFunction* K80RateMatrixFunction::clone( void ) const
{
    
    return new K80RateMatrixFunction( *this );
}


void K80RateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    double k = kappa->getValue();
    
    // set the base frequencies
    static_cast< RateMatrix_Kimura80* >(value)->setKappa( k );
    
    value->update();
}



void K80RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}


