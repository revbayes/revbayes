#include "revPoMoThree4NRateMatrixFunction.h"
#include "RateMatrix_revPoMoThree4N.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoThree4NRateMatrixFunction::revPoMoThree4NRateMatrixFunction( const TypedDagNode< std::int64_t > *n, 
                                                                    const TypedDagNode< Simplex > *bf,
                                                                    const TypedDagNode< RbVector<double> > *ex ,
                                                                    const TypedDagNode< RbVector<double> > *fc   ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoThree4N() ), 
N( n ), 
pi( bf), 
rho( ex ),
phi( fc )

{
   // add the lambda parameter as a parent
   addParameter( N );
   addParameter( pi );
   addParameter( rho );
   addParameter( phi );

   update();
}


revPoMoThree4NRateMatrixFunction::~revPoMoThree4NRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoThree4NRateMatrixFunction* revPoMoThree4NRateMatrixFunction::clone( void ) const
{
    return new revPoMoThree4NRateMatrixFunction( *this );
}


void revPoMoThree4NRateMatrixFunction::update( void )
{   
    double n = N ->getValue();
    const std::vector<double>& bf = pi->getValue();
    const std::vector<double>& ex = rho->getValue();
    const std::vector<double>& fc = phi->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoThree4N* >(value)->setN( n );
    static_cast< RateMatrix_revPoMoThree4N* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMoThree4N* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMoThree4N* >(value)->setPhi( fc );
    static_cast< RateMatrix_revPoMoThree4N* >(value)->update();
}


void revPoMoThree4NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N )
    {
        N = static_cast< const TypedDagNode< std::int64_t >* >( newP );
    }

    else if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    else if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}
