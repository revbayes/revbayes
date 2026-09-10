#include "PoMo4NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMo4NRateMatrixFunction::PoMo4NRateMatrixFunction(   const TypedDagNode< std::int64_t > *ni, 
                                                      const TypedDagNode< RbVector<double> > *m, 
                                                      const TypedDagNode< RbVector<double> > *f ) : 
TypedFunction<RateGenerator>( new RateMatrix_PoMo4N( computeNumStates( ni->getValue() ), ni->getValue() ) ),
N( ni ),
mu( m ),
phi( f )
{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( mu );
    addParameter( phi );

    update();
}


PoMo4NRateMatrixFunction::~PoMo4NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


std::int64_t PoMo4NRateMatrixFunction::computeNumStates( std::int64_t ni )
{

    std::int64_t numStates = 4 + 6*(ni-1);
    
    return numStates;

}




PoMo4NRateMatrixFunction* PoMo4NRateMatrixFunction::clone( void ) const
{
    
    return new PoMo4NRateMatrixFunction( *this );
}


void PoMo4NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    std::int64_t ni = N->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& f = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMo4N* >(value)->setN( ni );
    static_cast< RateMatrix_PoMo4N* >(value)->setMu( m );
    static_cast< RateMatrix_PoMo4N* >(value)->setPhi( f );

    value->update();
}



void PoMo4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
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

}




