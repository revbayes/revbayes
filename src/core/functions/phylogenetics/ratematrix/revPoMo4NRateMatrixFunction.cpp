#include "revPoMo4NRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMo4NRateMatrixFunction::revPoMo4NRateMatrixFunction(    const TypedDagNode< std::int64_t > *ni, 
                                                             const TypedDagNode< Simplex > *bf,
                                                             const TypedDagNode< RbVector<double> > *ex, 
                                                             const TypedDagNode< RbVector<double> > *f ) : 
TypedFunction<RateGenerator>( new RateMatrix_revPoMo4N( computeNumStates( ni->getValue() ),  ni->getValue()  ) ),
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


revPoMo4NRateMatrixFunction::~revPoMo4NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


std::int64_t revPoMo4NRateMatrixFunction::computeNumStates( std::int64_t ni )
{

    std::int64_t numStates = 4+6*(ni-1);
    
    return numStates;

}




revPoMo4NRateMatrixFunction* revPoMo4NRateMatrixFunction::clone( void ) const
{
    
    return new revPoMo4NRateMatrixFunction( *this );
}


void revPoMo4NRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    std::int64_t                        ni = N->getValue();
    const Simplex&              bf = pi->getValue();
    const std::vector<double>&  ex = rho->getValue();
    const std::vector<double>&  f  = phi->getValue();

    // set the base frequencies
    static_cast< RateMatrix_revPoMo4N* >(value)->setN( ni );
    static_cast< RateMatrix_revPoMo4N* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMo4N* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMo4N* >(value)->setPhi( f );

    value->update();
}



void revPoMo4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == N)
    {
        N =  static_cast<const TypedDagNode< std::int64_t >* >( newP );
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




