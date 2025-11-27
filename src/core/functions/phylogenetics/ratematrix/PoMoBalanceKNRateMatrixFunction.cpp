#include "PoMoBalanceKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoBalanceKNRateMatrixFunction::PoMoBalanceKNRateMatrixFunction(    const TypedDagNode< std::int64_t > *na,
                                                                     const TypedDagNode< std::int64_t > *ni,
                                                                     const TypedDagNode< RbVector<double> > *m,
                                                                     const TypedDagNode< RbVector<double> > *f,
                                                                     const TypedDagNode< RbVector<double> > *b,
                                                                     const TypedDagNode< RbVector<std::int64_t> > *Bf  ) :
TypedFunction<RateGenerator>( new RateMatrix_PoMoBalanceKN( computeNumStates( na->getValue(), ni->getValue() ), na->getValue(), ni->getValue(), computeNumMutRates( na->getValue() ) ) ),
K( na ),
N( ni ),
mu( m ),
phi( f ),
beta( b ),
B( Bf )

{
    // add the lambda parameter as a parent
    addParameter( K );
    addParameter( N );
    addParameter( mu );
    addParameter( phi );
    addParameter( beta );
    addParameter( B );

    update();
}


PoMoBalanceKNRateMatrixFunction::~PoMoBalanceKNRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


std::int64_t PoMoBalanceKNRateMatrixFunction::computeNumStates( std::int64_t na, std::int64_t ni )
{

    std::int64_t numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



std::int64_t PoMoBalanceKNRateMatrixFunction::computeNumMutRates( std::int64_t na )
{

    std::int64_t numMutRates = na*na-na;

    return numMutRates;

}


PoMoBalanceKNRateMatrixFunction* PoMoBalanceKNRateMatrixFunction::clone( void ) const
{
    
    return new PoMoBalanceKNRateMatrixFunction( *this );
}


void PoMoBalanceKNRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    std::int64_t na = K->getValue();
    std::int64_t ni = N->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& f = phi->getValue();
    const std::vector<double>& b  = beta->getValue();
    const std::vector<std::int64_t>&   Bf = B->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setK( na );
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setN( ni );
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setMu( m );
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setPhi( f );
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setBeta( b );
    static_cast< RateMatrix_PoMoBalanceKN* >(value)->setB( Bf );

    value->update();
}



void PoMoBalanceKNRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == K)
    {
        K = static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }
    
    if (oldP == N)
    {
        N =  static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }

    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
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




