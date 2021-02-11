#include "revPoMoKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoKNRateMatrixFunction::revPoMoKNRateMatrixFunction(   const TypedDagNode< long > *na, 
                                                      const TypedDagNode< long > *ni, 
                                                      const TypedDagNode< Simplex > *bf,
                                                      const TypedDagNode< RbVector<double> > *ex, 
                                                      const TypedDagNode< RbVector<double> > *f ) : 
TypedFunction<RateGenerator>( new RateMatrix_revPoMoKN( computeNumStates(na->getValue(), ni->getValue() ), na->getValue(), ni->getValue() , computeNumExchangeabilities( na->getValue() ) ) ),
K( na ),
N( ni),
pi( bf ),
rho( ex ),
phi( f )
{
    // add the lambda parameter as a parent
    addParameter( K );
    addParameter( N );
    addParameter( pi );
    addParameter( rho );
    addParameter( phi );

    update();
}


revPoMoKNRateMatrixFunction::~revPoMoKNRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long revPoMoKNRateMatrixFunction::computeNumStates( long na, long ni )
{

    long numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



long revPoMoKNRateMatrixFunction::computeNumExchangeabilities( long na )
{

    long numExchangeabilities = (na*na-na)*0.5;
    
    return numExchangeabilities;

}


revPoMoKNRateMatrixFunction* revPoMoKNRateMatrixFunction::clone( void ) const
{
    
    return new revPoMoKNRateMatrixFunction( *this );
}


void revPoMoKNRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    long                        na = K->getValue();
    long                        ni = N->getValue();
    const Simplex&              bf = pi->getValue();
    const std::vector<double>&  ex = rho->getValue();
    const std::vector<double>&  f  = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMoKN* >(value)->setK( na );
    static_cast< RateMatrix_revPoMoKN* >(value)->setN( ni );
    static_cast< RateMatrix_revPoMoKN* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMoKN* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMoKN* >(value)->setPhi( f );

    value->update();
}



void revPoMoKNRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == K)
    {
        K = static_cast<const TypedDagNode< long >* >( newP );
    }
    
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
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}




