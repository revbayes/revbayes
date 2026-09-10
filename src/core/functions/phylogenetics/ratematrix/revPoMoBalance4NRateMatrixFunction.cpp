#include "revPoMoBalance4NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoBalance4NRateMatrixFunction::revPoMoBalance4NRateMatrixFunction(   const TypedDagNode< std::int64_t > *n,
                                                                const TypedDagNode< Simplex  > *p,
                                                                const TypedDagNode< RbVector<double> > *r,
                                                                const TypedDagNode< RbVector<double> > *s,
                                                                const TypedDagNode< RbVector<double> > *b,
                                                                const TypedDagNode< RbVector<std::int64_t> > *Bf     ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoBalance4N( 4+6*(n->getValue()-1), n->getValue() ) ),
N( n ),
pi( p ),
rho( r ),
phi( s ),
beta( b ),
B( Bf )

{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( pi );
    addParameter( rho );
    addParameter( phi );
    addParameter( beta );
    addParameter( B );


    update();
}


revPoMoBalance4NRateMatrixFunction::~revPoMoBalance4NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoBalance4NRateMatrixFunction* revPoMoBalance4NRateMatrixFunction::clone( void ) const
{
    
    return new revPoMoBalance4NRateMatrixFunction( *this );
}


void revPoMoBalance4NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    std::int64_t n = N->getValue();
    const std::vector<double>& p  = pi->getValue();
    const std::vector<double>& r  = rho->getValue();
    const std::vector<double>& s  = phi->getValue();
    const std::vector<double>& b  = beta->getValue();
    const std::vector<std::int64_t>&   Bf = B->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setPi( p );
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setRho( r );
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setPhi( s );
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setBeta( b );
    static_cast< RateMatrix_revPoMoBalance4N* >(value)->setB( Bf );

    value->update();
}



void revPoMoBalance4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == N)
    { 
        N = static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }
    
    if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == B)
    {
        B = static_cast<const TypedDagNode< RbVector<std::int64_t> >* >( newP );
    }


}




