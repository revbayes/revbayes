//
//  CovarionRateMatrixFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "CovarionRateMatrixFunction.h"

#include <cstddef>

#include "RateMatrix_Covarion.h"
#include "RateGenerator.h"
#include "TypedFunction.h"
#include "AbstractRateMatrix.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

CovarionRateMatrixFunction::CovarionRateMatrixFunction(const TypedDagNode< RbVector<RateGenerator> > *rm,
                                                       const TypedDagNode< RateGenerator > *sr,
                                                       const TypedDagNode< RbVector<double> > *cr,
                                                       bool rescaled) :
TypedFunction<RateGenerator>( new RateMatrix_Covarion( rm->getValue()[0].size(), cr->getValue().size(), rescaled ) ),
rate_matrices( rm ),
clock_rates( cr ),
switch_rates( sr )

{
    // add the lambda parameter as a parent
    addParameter( rate_matrices );
    addParameter( clock_rates );
    addParameter( switch_rates );
    
    update();
}

CovarionRateMatrixFunction::~CovarionRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CovarionRateMatrixFunction* CovarionRateMatrixFunction::clone( void ) const
{
    return new CovarionRateMatrixFunction( *this );
}


void CovarionRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    //    const RbVector<RateMatrix>& rm = rate_matrices->getValue();
    RbVector<MatrixReal> rm;
    for (size_t i = 0; i < rate_matrices->getValue().size(); i++)
    {
        const AbstractRateMatrix* rm_ptr = dynamic_cast<const AbstractRateMatrix*>( &rate_matrices->getValue()[i] );
        rm.push_back( rm_ptr->getRateMatrix() );
    }
    const AbstractRateMatrix* sr_ptr = dynamic_cast<const AbstractRateMatrix*>( &switch_rates->getValue() );
    const MatrixReal& sr           = sr_ptr->getRateMatrix();
    const RbVector<double>& cr     = clock_rates->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_Covarion* >(value)->setRateMatrices( rm );
    static_cast< RateMatrix_Covarion* >(value)->setSwitchRates( sr );
    static_cast< RateMatrix_Covarion* >(value)->setClockRates( cr );
    
    value->update();
    
}



void CovarionRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == rate_matrices)
    {
        rate_matrices = static_cast<const TypedDagNode< RbVector<RateGenerator> >* >( newP );
    }
    
    if (oldP == switch_rates)
    {
        switch_rates = static_cast<const TypedDagNode<RateGenerator>* >( newP );
    }
    
    if (oldP == clock_rates)
    {
        clock_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

