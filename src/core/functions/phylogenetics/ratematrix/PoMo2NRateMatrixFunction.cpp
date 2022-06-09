#include "PoMo2NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMo2NRateMatrixFunction::PoMo2NRateMatrixFunction(   long virtual_N,
                                                      const TypedDagNode< double > *ni,
                                                      const TypedDagNode< RbVector<double> > *m, 
                                                      const TypedDagNode< RbVector<double> > *f,
                                                      bool use_mutation_correction,
                                                      bool use_drift_correction) :
TypedFunction<RateGenerator>( new RateMatrix_PoMo2N( computeNumStates( virtual_N ), virtual_N, use_mutation_correction, use_drift_correction ) ),
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


PoMo2NRateMatrixFunction::~PoMo2NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long PoMo2NRateMatrixFunction::computeNumStates( long ni )
{

    long numStates = ni+1;
    
    return numStates;

}




PoMo2NRateMatrixFunction* PoMo2NRateMatrixFunction::clone( void ) const
{
    
    return new PoMo2NRateMatrixFunction( *this );
}


void PoMo2NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    double ni = N->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& f = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_PoMo2N* >(value)->setNeff( ni );
    static_cast< RateMatrix_PoMo2N* >(value)->setMu( m );
    static_cast< RateMatrix_PoMo2N* >(value)->setPhi( f );

    value->update();
}



void PoMo2NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == N)
    {
        N =  static_cast<const TypedDagNode< double >* >( newP );
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




