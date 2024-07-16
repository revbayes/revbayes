#include "PoMoKNrecurrentMutationsRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoKNrecurrentMutationsRateMatrixFunction::PoMoKNrecurrentMutationsRateMatrixFunction(     long na,
                                                                                            long nv,
                                                                                            const TypedDagNode< RbVector<double> > *mut,
                                                                                            const TypedDagNode< RbVector<double> > *fit,
                                                                                            bool rm ) :
TypedFunction<RateGenerator>( new RateMatrix_PoMoKNrecurrentMutations( computeNumStates(na, nv), na, nv, computeNumMutRates( na ), rm ) ),
K( na ),
V( nv ),
mu( mut ),
phi( fit ),
R( rm )
{
    // add the lambda parameter as a parent
    addParameter( mu );
    addParameter( phi );

    update();
}


PoMoKNrecurrentMutationsRateMatrixFunction::~PoMoKNrecurrentMutationsRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


long PoMoKNrecurrentMutationsRateMatrixFunction::computeNumStates( long na, long ni )
{

    long numStates = na + (na*na-na)*(ni-1)*0.5;
    
    return numStates;

}



long PoMoKNrecurrentMutationsRateMatrixFunction::computeNumMutRates( long na )
{

    long numMutRates = na*na-na;
    
    return numMutRates;

}


PoMoKNrecurrentMutationsRateMatrixFunction* PoMoKNrecurrentMutationsRateMatrixFunction::clone( void ) const
{
    
    return new PoMoKNrecurrentMutationsRateMatrixFunction( *this );
}


void PoMoKNrecurrentMutationsRateMatrixFunction::update( void )
{

    if ( mu != NULL )
    {
        static_cast< RateMatrix_PoMoKNrecurrentMutations* >(value)->setMu( mu->getValue() );
    }
    if ( phi != NULL )
    {
        static_cast< RateMatrix_PoMoKNrecurrentMutations* >(value)->setPhi( phi->getValue() );
    }
    
    value->update();
}



void PoMoKNrecurrentMutationsRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == phi)
    {
        phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}




