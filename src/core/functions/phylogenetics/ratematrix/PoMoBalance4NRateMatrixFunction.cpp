#include "PoMoBalance4NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoBalance4NRateMatrixFunction::PoMoBalance4NRateMatrixFunction(   const TypedDagNode< long > *ni,
                                                      const TypedDagNode< RbVector<double> > *m, 
                                                      const TypedDagNode< RbVector<double> > *f,
                                                      const TypedDagNode< RbVector<double> > *b,
                                                      const TypedDagNode< RbVector<long> > *Bf  ) :
TypedFunction<RateGenerator>( new RateMatrix_PoMoBalance4N( computeNumStates( ni->getValue() ), ni->getValue() ) ),
N( ni ),
mu( m ),
phi( f ),
beta( b ),
B( Bf )

{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( mu );
    addParameter( phi );
    addParameter( beta );
    addParameter( B );

    update();
}


PoMoBalance4NRateMatrixFunction::~PoMoBalance4NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long PoMoBalance4NRateMatrixFunction::computeNumStates( long ni )
{

    long numStates = 4 + 6*(ni-1);
    
    return numStates;

}




PoMoBalance4NRateMatrixFunction* PoMoBalance4NRateMatrixFunction::clone( void ) const
{
    
    return new PoMoBalance4NRateMatrixFunction( *this );
}


void PoMoBalance4NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    long ni = N->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& f = phi->getValue();
    const std::vector<double>& b  = beta->getValue();
    const std::vector<long>&   Bf = B->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMoBalance4N* >(value)->setN( ni );
    static_cast< RateMatrix_PoMoBalance4N* >(value)->setMu( m );
    static_cast< RateMatrix_PoMoBalance4N* >(value)->setPhi( f );
    static_cast< RateMatrix_PoMoBalance4N* >(value)->setBeta( b );
    static_cast< RateMatrix_PoMoBalance4N* >(value)->setB( Bf );

    value->update();
}



void PoMoBalance4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
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

    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == B)
    {
        B = static_cast<const TypedDagNode< RbVector<long> >* >( newP );
    }

}




