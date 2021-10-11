#include "revPoMo2NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMo2NRateMatrixFunction::revPoMo2NRateMatrixFunction(    const TypedDagNode< long > *ni, 
                                                             const TypedDagNode< Simplex > *bf,
                                                             const TypedDagNode< double > *ex, 
                                                             const TypedDagNode< RbVector<double> > *f ) : 
TypedFunction<RateGenerator>( new RateMatrix_revPoMo2N( computeNumStates( ni->getValue() ),  ni->getValue()  ) ),
N( ni),
pi( bf ),
rho( ex ),
phi( f )
{
    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( pi );
    addParameter( rho );
    addParameter( phi );

    update();
}


revPoMo2NRateMatrixFunction::~revPoMo2NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long revPoMo2NRateMatrixFunction::computeNumStates( long ni )
{

    long numStates = ni+1;
    
    return numStates;

}




revPoMo2NRateMatrixFunction* revPoMo2NRateMatrixFunction::clone( void ) const
{
    
    return new revPoMo2NRateMatrixFunction( *this );
}


void revPoMo2NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    long                        ni = N->getValue();
    const Simplex&              bf = pi->getValue();
    double                      ex = rho->getValue();
    const std::vector<double>&  f  = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMo2N* >(value)->setN( ni );
    static_cast< RateMatrix_revPoMo2N* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMo2N* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMo2N* >(value)->setPhi( f );

    value->update();
}



void revPoMo2NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

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
        rho = static_cast<const TypedDagNode< double >* >( newP );
    }

    if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}




