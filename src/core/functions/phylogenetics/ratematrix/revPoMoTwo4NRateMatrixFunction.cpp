#include "revPoMoTwo4NRateMatrixFunction.h"
#include "RateMatrix_revPoMoTwo4N.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoTwo4NRateMatrixFunction::revPoMoTwo4NRateMatrixFunction( const TypedDagNode< std::int64_t > *n, 
                                                                const TypedDagNode< Simplex > *bf,
                                                                const TypedDagNode< RbVector<double> > *ex  ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoTwo4N() ), 
N( n ), 
pi( bf), 
rho( ex )

{
   // add the lambda parameter as a parent
   addParameter( N );
   addParameter( pi );
   addParameter( rho );

   update();
}


revPoMoTwo4NRateMatrixFunction::~revPoMoTwo4NRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoTwo4NRateMatrixFunction* revPoMoTwo4NRateMatrixFunction::clone( void ) const
{
    return new revPoMoTwo4NRateMatrixFunction( *this );
}


void revPoMoTwo4NRateMatrixFunction::update( void )
{   
    double n = N ->getValue();
    const std::vector<double>& bf = pi->getValue();
    const std::vector<double>& ex = rho->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoTwo4N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoTwo4N* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMoTwo4N* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMoTwo4N* >(value)->update();
}


void revPoMoTwo4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N )
    {
        N = static_cast< const TypedDagNode< std::int64_t >* >( newP );
    }

    else if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }


}
