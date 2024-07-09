#include "revPoMo2NrecurrentRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMo2NrecurrentRateMatrixFunction::revPoMo2NrecurrentRateMatrixFunction(  const TypedDagNode< long > *ni, 
                                                      const TypedDagNode< Simplex > *bf,
                                                      const TypedDagNode< double > *ex ) : 
TypedFunction<RateGenerator>( new RateMatrix_revPoMo2Nrecurrent( ni->getValue() + 1, ni->getValue() ) ),
N( ni),
pi( bf ),
rho( ex )
{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( pi );
    addParameter( rho );

    update();
}


revPoMo2NrecurrentRateMatrixFunction::~revPoMo2NrecurrentRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



revPoMo2NrecurrentRateMatrixFunction* revPoMo2NrecurrentRateMatrixFunction::clone( void ) const
{
    
    return new revPoMo2NrecurrentRateMatrixFunction( *this );
}


void revPoMo2NrecurrentRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    long                        ni = N->getValue();
    const Simplex&              bf = pi->getValue();
    double                      ex = rho->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMo2Nrecurrent* >(value)->setN( ni );
    static_cast< RateMatrix_revPoMo2Nrecurrent* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMo2Nrecurrent* >(value)->setRho( ex );

    value->update();
}



void revPoMo2NrecurrentRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == N)
    {
        N =  static_cast<const TypedDagNode< long >* >( newP );
    }

    if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< double >* >( newP );
    }


}




