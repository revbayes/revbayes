#include "revPoMoNeutralM4NRateMatrixFunction.h"
#include "RateMatrix_revPoMoNeutralM4N.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoNeutralM4NRateMatrixFunction::revPoMoNeutralM4NRateMatrixFunction( const TypedDagNode< std::int64_t > *n, 
                                                                          const TypedDagNode< std::int64_t > *m, 
                                                                          const TypedDagNode< Simplex > *bf,
                                                                          const TypedDagNode< RbVector<double> > *ex  ) : 
TypedFunction<RateGenerator>( new RateMatrix_revPoMoNeutralM4N( computeNumStates( m->getValue() ),  m->getValue()  ) ), 
N( n ),
M( m ), 
pi( bf), 
rho( ex )

{
   // add the lambda parameter as a parent
   addParameter( N );
   addParameter( M );
   addParameter( pi );
   addParameter( rho );

   update();
}


revPoMoNeutralM4NRateMatrixFunction::~revPoMoNeutralM4NRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


std::int64_t revPoMoNeutralM4NRateMatrixFunction::computeNumStates( std::int64_t mi )
{

    std::int64_t numStates = 4+6*(mi-1);
    
    return numStates;

}


revPoMoNeutralM4NRateMatrixFunction* revPoMoNeutralM4NRateMatrixFunction::clone( void ) const
{
    return new revPoMoNeutralM4NRateMatrixFunction( *this );
}


void revPoMoNeutralM4NRateMatrixFunction::update( void )
{   
    std::int64_t n = N ->getValue();
    std::int64_t m = M ->getValue();  
    const std::vector<double>& bf = pi->getValue();
    const std::vector<double>& ex = rho->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoNeutralM4N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoNeutralM4N* >(value)->setM( m );
    static_cast< RateMatrix_revPoMoNeutralM4N* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMoNeutralM4N* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMoNeutralM4N* >(value)->update();
}


void revPoMoNeutralM4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N )
    {
        N = static_cast< const TypedDagNode< std::int64_t >* >( newP );
    }

    else if (oldP == M )
    {
        M = static_cast< const TypedDagNode< std::int64_t >* >( newP );
    }

    else if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }


}
