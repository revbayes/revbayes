#include "PoMoKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoKNRateMatrixFunction::PoMoKNRateMatrixFunction(   long num_all,
                                                      long virt,
                                                      const TypedDagNode< double > *eff,
                                                      const TypedDagNode< RbVector<double> > *mut,
                                                      const TypedDagNode< RbVector<double> > *sel ) :
TypedFunction<RateGenerator>( new RateMatrix_PoMoKN( computeNumStates(num_all, virt), num_all, virt, virt, computeNumMutRates( num_all ) ) ),
num_alleles( num_all ),
virt_pop_size( virt ),
eff_pop_size( eff ),
mu( mut ),
phi( sel )
{
    // add the lambda parameter as a parent
    addParameter( eff_pop_size );
    addParameter( mu );
    addParameter( phi );

    update();
}


PoMoKNRateMatrixFunction::~PoMoKNRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long PoMoKNRateMatrixFunction::computeNumStates( long na, long ni )
{

    long numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



long PoMoKNRateMatrixFunction::computeNumMutRates( long na )
{

    long numMutRates = na*na-na;
    
    return numMutRates;

}


PoMoKNRateMatrixFunction* PoMoKNRateMatrixFunction::clone( void ) const
{
    
    return new PoMoKNRateMatrixFunction( *this );
}


void PoMoKNRateMatrixFunction::update( void )
{

    if ( eff_pop_size != NULL )
    {
        static_cast< RateMatrix_PoMoKN* >(value)->setEffectivePopulationSize( eff_pop_size->getValue() );
    }
    if ( mu != NULL )
    {
        static_cast< RateMatrix_PoMoKN* >(value)->setMu( mu->getValue() );
    }
    if ( phi != NULL )
    {
        static_cast< RateMatrix_PoMoKN* >(value)->setPhi( phi->getValue() );
    }
    
    value->update();
}



void PoMoKNRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == eff_pop_size)
    {
        eff_pop_size =  static_cast<const TypedDagNode< double >* >( newP );
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




