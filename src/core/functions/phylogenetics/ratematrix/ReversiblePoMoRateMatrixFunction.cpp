#include "ReversiblePoMoRateMatrixFunction.h"
#include "RateMatrix_PoMo.h"
#include "RbException.h"

using namespace RevBayesCore;


//  //MJL 140822: caused compile error
// ReversiblePoMoRateMatrixFunction::ReversiblePoMoRateMatrixFunction(const TypedDagNode< int > *ps, const TypedDagNode< RateGenerator > *mm ) : TypedFunction<RateGenerator>( new RateMatrix_ReversiblePoMo(4 + 6*(ps->getValue() - 1), mm->getValue(), ps->getValue() ) ), populationSize( ps ), mutationMatrix( mm ) {
//     // add the lambda parameter as a parent
//     addParameter( populationSize );
//     addParameter( mutationMatrix );
//
//     update();
// }


ReversiblePoMoRateMatrixFunction::ReversiblePoMoRateMatrixFunction(const TypedDagNode< long > *ps, const TypedDagNode< RbVector<double> > *rho, const TypedDagNode< Simplex > *pi   ) : TypedFunction<RateGenerator>( new RateMatrix_ReversiblePoMo(4 + 6*(ps->getValue() - 1), rho->getValue(), pi->getValue(), ps->getValue() ) ), populationSize( ps ), exchangeabilities( rho ), equilibriumFrequencies( pi ) {
   // add the lambda parameter as a parent
   addParameter( populationSize );
   addParameter( exchangeabilities );
   addParameter( equilibriumFrequencies );

   update();
}



ReversiblePoMoRateMatrixFunction::~ReversiblePoMoRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ReversiblePoMoRateMatrixFunction* ReversiblePoMoRateMatrixFunction::clone( void ) const
{
    return new ReversiblePoMoRateMatrixFunction( *this );
}


void ReversiblePoMoRateMatrixFunction::update( void )
{

    // set the mutation rates
    //static_cast< RateMatrix_ReversiblePoMo* >(value)->setMutationRates( static_cast< const RateMatrix& > ( mutationMatrix->getValue() ) );
    static_cast< RateMatrix_ReversiblePoMo* >(value)->update();

}



void ReversiblePoMoRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    // if (oldP == mutationMatrix)
    // {
    //     mutationMatrix = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    // }
    if (oldP == exchangeabilities)
    {
        exchangeabilities = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == equilibriumFrequencies)
    {
        equilibriumFrequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

}
