#include "PoMoBalanceRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoBalanceRateMatrixFunction::PoMoBalanceRateMatrixFunction(   const TypedDagNode< double > *n,
                                                                const TypedDagNode< Simplex  > *p,
                                                                const TypedDagNode< RbVector<double> > *r,
                                                                const TypedDagNode< RbVector<double> > *s,
                                                                const TypedDagNode< double > *b ) : TypedFunction<RateGenerator>( new RateMatrix_PoMoBalance( 4+6*(n->getValue()-1), n->getValue() ) ),
N( n ),
pi( p ),
rho( r ),
sigma( s ),
beta( b )

{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( pi );
    addParameter( rho );
    addParameter( sigma );
    addParameter( beta );


    update();
}


PoMoBalanceRateMatrixFunction::~PoMoBalanceRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PoMoBalanceRateMatrixFunction* PoMoBalanceRateMatrixFunction::clone( void ) const
{
    
    return new PoMoBalanceRateMatrixFunction( *this );
}


void PoMoBalanceRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    double n = N->getValue();
    const std::vector<double>& p = pi->getValue();
    const std::vector<double>& r = rho->getValue();
    const std::vector<double>& s = sigma->getValue();
    double b = beta->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMoBalance* >(value)->setN( n );
    static_cast< RateMatrix_PoMoBalance* >(value)->setPi( p );
    static_cast< RateMatrix_PoMoBalance* >(value)->setRho( r );
    static_cast< RateMatrix_PoMoBalance* >(value)->setSigma( s );
    static_cast< RateMatrix_PoMoBalance* >(value)->setBeta( b );

    value->update();
}



void PoMoBalanceRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == N)
    {
        N = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode< double >* >( newP );
    }

}




