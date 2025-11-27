#include "PoMoKNRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoKNRateMatrixFunction::PoMoKNRateMatrixFunction(   std::int64_t num_all,
                                                      std::int64_t virt,
                                                      const TypedDagNode< double > *eff,
                                                      const TypedDagNode< RbVector<double> > *mut,
                                                      const TypedDagNode< RbVector<double> > *sel ) :
TypedFunction<RateGenerator>( new RateMatrix_PoMoKN( computeNumStates(num_all, virt), num_all, virt, (eff != NULL ? eff->getValue() : virt), computeNumMutRates( num_all ) ) ),
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


std::int64_t PoMoKNRateMatrixFunction::computeNumStates( std::int64_t na, std::int64_t ni )
{

    std::int64_t numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



std::int64_t PoMoKNRateMatrixFunction::computeNumMutRates( std::int64_t na )
{

    std::int64_t numMutRates = na*na-na;
    
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
