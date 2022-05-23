#include "revPoMoM2NRateMatrixFunction.h"
#include "RateMatrix_revPoMoM2N.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoM2NRateMatrixFunction::revPoMoM2NRateMatrixFunction( const TypedDagNode< long > *n, 
                                                            const TypedDagNode< RbVector<double> > *m ,
                                                            const TypedDagNode< long > *v, 
                                                            const TypedDagNode< double > *g ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoM2N( v->getValue() + 1 ) ), 
N( n ), 
mu( m ),
M( v ),
gen( g )
{
   // add the lambda parameter as a parent
   addParameter( N );
   addParameter( mu );
   addParameter( M );
   addParameter( gen );

   update();
}


revPoMoM2NRateMatrixFunction::~revPoMoM2NRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoM2NRateMatrixFunction* revPoMoM2NRateMatrixFunction::clone( void ) const
{
    return new revPoMoM2NRateMatrixFunction( *this );
}


void revPoMoM2NRateMatrixFunction::update( void )
{   
    long n   = N ->getValue();
    const std::vector<double>& m = mu->getValue();
    long v   = M ->getValue();
    double g = gen ->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoM2N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoM2N* >(value)->setMu( m );
    static_cast< RateMatrix_revPoMoM2N* >(value)->setM( v );
    static_cast< RateMatrix_revPoMoM2N* >(value)->setGen( g );
    static_cast< RateMatrix_revPoMoM2N* >(value)->update();
}


void revPoMoM2NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N )
    {
        N = static_cast< const TypedDagNode< long >* >( newP );
    }

    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    else if (oldP == M)
    {
        M = static_cast<const TypedDagNode< long >* >( newP );
    }

    else if (oldP == gen)
    {
        gen = static_cast<const TypedDagNode< double >* >( newP );
    }



}
