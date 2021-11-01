#include "revPoMoThree4RateMatrixFunction.h"
#include "RateMatrix_revPoMoThree4.h"
#include "RbException.h"

using namespace RevBayesCore;


revPoMoThree4RateMatrixFunction::revPoMoThree4RateMatrixFunction(   const TypedDagNode< Simplex > *bf,
                                                                    const TypedDagNode< RbVector<double> > *ex ,
                                                                    const TypedDagNode< RbVector<double> > *fc   ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoThree4() ), 
pi( bf), 
rho( ex ),
phi( fc )

{
   // add the lambda parameter as a parent
   addParameter( pi );
   addParameter( rho );
   addParameter( phi );

   update();
}


revPoMoThree4RateMatrixFunction::~revPoMoThree4RateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoThree4RateMatrixFunction* revPoMoThree4RateMatrixFunction::clone( void ) const
{
    return new revPoMoThree4RateMatrixFunction( *this );
}


void revPoMoThree4RateMatrixFunction::update( void )
{   
    const std::vector<double>& bf = pi->getValue();
    const std::vector<double>& ex = rho->getValue();
    const std::vector<double>& fc = phi->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoThree4* >(value)->setPi( bf );
    static_cast< RateMatrix_revPoMoThree4* >(value)->setRho( ex );
    static_cast< RateMatrix_revPoMoThree4* >(value)->setPhi( fc );
    static_cast< RateMatrix_revPoMoThree4* >(value)->update();
}


void revPoMoThree4RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
   
    if (oldP == pi)
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
