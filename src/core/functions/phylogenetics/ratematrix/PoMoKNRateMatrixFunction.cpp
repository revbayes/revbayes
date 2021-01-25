#include "PoMoKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoKNRateMatrixFunction::PoMoKNRateMatrixFunction(   const TypedDagNode< long > *na, 
                                                      const TypedDagNode< long > *ni, 
                                                      const TypedDagNode< RbVector<double> > *m, 
                                                      const TypedDagNode< RbVector<double> > *f ) : 
TypedFunction<RateGenerator>( new RateMatrix_PoMoKN( computeNumStates(na->getValue(), ni->getValue() ), na->getValue(), ni->getValue() , computeNumMutRates( na->getValue() ) ) ),
K( na ),
N( ni),
mu( m ),
phi( f )
{
    // add the lambda parameter as a parent
    addParameter( K );
    addParameter( N );
    addParameter( mu );
    addParameter( phi );

    update();
}


PoMoKNRateMatrixFunction::~PoMoKNRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long PoMoKNRateMatrixFunction::computeNumStates( long na, long ni )
{

    long numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



long PoMoKNRateMatrixFunction::computeNumMutRates( long na )
{

    long numMutRates = na*na-na;
    
    return numMutRates;

}


PoMoKNRateMatrixFunction* PoMoKNRateMatrixFunction::clone( void ) const
{
    
    return new PoMoKNRateMatrixFunction( *this );
}


void PoMoKNRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    long na = K->getValue();
    long ni = N->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& f = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMoKN* >(value)->setK( na );
    static_cast< RateMatrix_PoMoKN* >(value)->setN( ni );
    static_cast< RateMatrix_PoMoKN* >(value)->setMu( m );
    static_cast< RateMatrix_PoMoKN* >(value)->setPhi( f );

    value->update();
}



void PoMoKNRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == K)
    {
        K = static_cast<const TypedDagNode< long >* >( newP );
    }
    
    if (oldP == N)
    {
        N =  static_cast<const TypedDagNode< long >* >( newP );
    }

    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}




