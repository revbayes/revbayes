#include "revPoMoM2NRateMatrixFunction.h"
#include "RateMatrix_revPoMoM2N.h"
#include "RbException.h"

using namespace RevBayesCore;

revPoMoM2NRateMatrixFunction::revPoMoM2NRateMatrixFunction( long m,
                                                            const TypedDagNode< double > *n,
                                                            const TypedDagNode< RbVector<double> > *mut ,
                                                            const TypedDagNode< double > *g ) : TypedFunction<RateGenerator>( new RateMatrix_revPoMoM2N( m ) ),
N_eff( n ),
mu( mut ),
N_virt( m ),
gen( g )
{
   // add the lambda parameter as a parent
   addParameter( N_eff );
   addParameter( mu );
   addParameter( gen );

   update();
}


revPoMoM2NRateMatrixFunction::~revPoMoM2NRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


revPoMoM2NRateMatrixFunction* revPoMoM2NRateMatrixFunction::clone( void ) const
{
    return new revPoMoM2NRateMatrixFunction( *this );
}


void revPoMoM2NRateMatrixFunction::update( void )
{   
    double n = N_eff->getValue();
    const std::vector<double>& m = mu->getValue();
    double g = gen->getValue();

    // set parameters
    static_cast< RateMatrix_revPoMoM2N* >(value)->setNEffective( n );
    static_cast< RateMatrix_revPoMoM2N* >(value)->setMu( m );
    static_cast< RateMatrix_revPoMoM2N* >(value)->setGen( g );
    static_cast< RateMatrix_revPoMoM2N* >(value)->update();
    
//    std::vector<double> phi = std::vector<double>(2,1.0);
//    static_cast< RateMatrix_PoMo2N* >(value)->setNeff( n );
//    static_cast< RateMatrix_PoMo2N* >(value)->setMu( m );
//    static_cast< RateMatrix_PoMo2N* >(value)->setPhi( phi );
//    static_cast< RateMatrix_PoMo2N* >(value)->update();
}


void revPoMoM2NRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == N_eff )
    {
        N_eff = static_cast< const TypedDagNode< double >* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == gen)
    {
        gen = static_cast<const TypedDagNode< double >* >( newP );
    }



}
