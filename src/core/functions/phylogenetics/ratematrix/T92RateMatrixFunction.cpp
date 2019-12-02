#include "T92RateMatrixFunction.h"
#include "RateMatrix_Tamura92.h"
#include "RbException.h"

using namespace RevBayesCore;

T92RateMatrixFunction::T92RateMatrixFunction(const TypedDagNode< double > *eqGc, const TypedDagNode< double > *tstv) : TypedFunction<RateGenerator>( new RateMatrix_Tamura92( ) ),
    equilibriumGc( eqGc ),
    transitionTransversionRate( tstv )
{

    // add the lambda parameter as a parent
    addParameter( equilibriumGc );
    addParameter( transitionTransversionRate );
    
    update();
}


T92RateMatrixFunction::~T92RateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



T92RateMatrixFunction* T92RateMatrixFunction::clone( void ) const
{
    return new T92RateMatrixFunction( *this );
}


void T92RateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    double k  = transitionTransversionRate->getValue();
    double gc = equilibriumGc->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_Tamura92* >(value)->setGC( gc );
    static_cast< RateMatrix_Tamura92* >(value)->setKappa( k );
    
    value->update();

}



void T92RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == equilibriumGc)
    {
        equilibriumGc = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == transitionTransversionRate)
    {
        transitionTransversionRate = static_cast<const TypedDagNode< double >* >( newP );
    }
}


