#include "revPoMoBalanceKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoBalanceKNRateMatrixFunction::revPoMoBalanceKNRateMatrixFunction(    const TypedDagNode<std::int64_t> *na,
                                                                const TypedDagNode<std::int64_t> *ni,
                                                                const TypedDagNode< Simplex  > *p,
                                                                const TypedDagNode< RbVector<double> > *r,
                                                                const TypedDagNode< RbVector<double> > *s,
                                                                const TypedDagNode< RbVector<double> > *b     ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoBalanceKN( computeNumStates(na->getValue(), ni->getValue() ), na->getValue(), ni->getValue() , computeNumExchangeabilities( na->getValue() ) ) ),

K( na ),
N( ni ),
pi( p ),
rho( r ),
phi( s ),
beta( b )

{
    // add the lambda parameter as a parent
    addParameter( K );
    addParameter( N );
    addParameter( pi );
    addParameter( rho );
    addParameter( phi );
    addParameter( beta );

    update();
}


revPoMoBalanceKNRateMatrixFunction::~revPoMoBalanceKNRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


std::int64_t revPoMoBalanceKNRateMatrixFunction::computeNumStates(std::int64_t na,std::int64_t ni )
{

    std::int64_t numStates = na + (na*na-na)*(ni-1)*0.5;

    return numStates;

}


std::int64_t revPoMoBalanceKNRateMatrixFunction::computeNumExchangeabilities(std::int64_t na )
{

    std::int64_t numExchangeabilities = (na*na-na)*0.5;

    return numExchangeabilities;

}


revPoMoBalanceKNRateMatrixFunction* revPoMoBalanceKNRateMatrixFunction::clone( void ) const
{
    
    return new revPoMoBalanceKNRateMatrixFunction( *this );
}


void revPoMoBalanceKNRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    std::int64_t                na = K->getValue();
    std::int64_t                ni = N->getValue();
    const std::vector<double>&  p  = pi->getValue();
    const std::vector<double>&  r  = rho->getValue();
    const std::vector<double>&  s  = phi->getValue();
    const std::vector<double>&  b  = beta->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setK( na );
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setN( ni );
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setPi( p );
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setRho( r );
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setPhi( s );
    static_cast< RateMatrix_revPoMoBalanceKN* >(value)->setBeta( b );

    value->update();
}



void revPoMoBalanceKNRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == K)
    {
        K = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
    if (oldP == N)
    { 
        N = static_cast<const TypedDagNode<std::int64_t>* >( newP );
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


}




