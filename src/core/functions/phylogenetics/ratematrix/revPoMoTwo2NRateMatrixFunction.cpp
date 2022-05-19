#include "revPoMoTwo2NRateMatrixFunction.h"
#include "RateMatrix_revPoMoTwo2N.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoTwo2NRateMatrixFunction::revPoMoTwo2NRateMatrixFunction( const TypedDagNode< long > *n, 
                                                                const TypedDagNode< RbVector<double> > *m  ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoTwo2N() ), 
N( n ), 
mu( m )

{
   // add the lambda parameter as a parent
   addParameter( N );
   addParameter( mu );

   update();
}


revPoMoTwo2NRateMatrixFunction::~revPoMoTwo2NRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoTwo2NRateMatrixFunction* revPoMoTwo2NRateMatrixFunction::clone( void ) const
{
    return new revPoMoTwo2NRateMatrixFunction( *this );
}


void revPoMoTwo2NRateMatrixFunction::update( void )
{   
    double n = N ->getValue();
    const std::vector<double>& m = mu->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoTwo2N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoTwo2N* >(value)->setMu( m );
    static_cast< RateMatrix_revPoMoTwo2N* >(value)->update();
}


void revPoMoTwo2NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N )
    {
        N = static_cast< const TypedDagNode< long >* >( newP );
    }

    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }


}
