//
//  CompositeRateMatrixFunction.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/13/24.
//

#include "CompositeRateMatrixFunction.h"

#include <cstddef>

//#include "RateMatrix_Composite.h"
//#include "RateMatrix_Composite.h"
#include "RateMatrix_FreeK.h"
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

CompositeRateMatrixFunction::CompositeRateMatrixFunction(const TypedDagNode< RbVector<RateGenerator> > *rm1, const TypedDagNode< RbVector<RateGenerator> > *rm2) :
TypedFunction<RateGenerator>( new RateMatrix_FreeK( rm1->getValue().size() * rm2->getValue().size() ) ),
rate_matrices1( rm1 ),
rate_matrices2( rm2 )

{
    // add the lambda parameter as a parent
    addParameter( rate_matrices1 );
    addParameter( rate_matrices2 );
//    addParameter( clock_rates );
//    addParameter( switch_rates );
    
    update();
}

CompositeRateMatrixFunction::~CompositeRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CompositeRateMatrixFunction* CompositeRateMatrixFunction::clone( void ) const
{
    return new CompositeRateMatrixFunction( *this );
}


/*
void CompositeRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    //    const RbVector<RateMatrix>& rm = rate_matrices->getValue();
    RbVector<MatrixReal> rm1;
    for (size_t i = 0; i < rate_matrices1->getValue().size(); i++)
    {
        const AbstractRateMatrix* rm_ptr = dynamic_cast<const AbstractRateMatrix*>( &rate_matrice1s->getValue()[i] );
        rm1.push_back( rm_ptr->getRateMatrix() );
    }
    RbVector<MatrixReal> rm2;
    for (size_t i = 0; i < rate_matrices2->getValue().size(); i++)
    {
        const AbstractRateMatrix* rm_ptr = dynamic_cast<const AbstractRateMatrix*>( &rate_matrices2->getValue()[i] );
        rm2.push_back( rm_ptr->getRateMatrix() );
    }
//    const AbstractRateMatrix* sr_ptr = dynamic_cast<const AbstractRateMatrix*>( &switch_rates->getValue() );
//    const MatrixReal& sr           = sr_ptr->getRateMatrix();
//    const RbVector<double>& cr     = clock_rates->getValue();
    
    
    // set the rate matrix values
    static_cast< RateMatrix_Composite* >(value)->setRateMatrices1( rm1 );
    static_cast< RateMatrix_Composite* >(value)->setRateMatrices2( rm2 );
//    static_cast< RateMatrix_Composite* >(value)->setSwitchRates( sr );
//    static_cast< RateMatrix_Composite* >(value)->setClockRates( cr );
    
    value->update();
    
}
*/

void CompositeRateMatrixFunction::update( void )
{
    size_t num_states_1 = rate_matrices2->getValue().size();    // number ranges == number trait mtx
    size_t num_states_2 = rate_matrices1->getValue().size();    // number traits == number biogeo mtx
    size_t num_states = num_states_1 * num_states_2;
    
    // set up a 2-d matrix to hold the rates for the combined observed and hidden states
    std::vector< std::vector<double> > rate_matrix = std::vector< std::vector<double> >( num_states, std::vector<double>( num_states, 0.0 ) );
    
    // populate rate matrix elements for rm1
    // diagonal blocks of arbitrary rate matrices
    for (size_t i = 0; i < num_states_2; i++) {
        size_t offset = i * num_states_1;
        for (size_t j = 0; j < num_states_1; j++) {
            for (size_t k = 0; k < num_states_1; k++) {
                rate_matrix[offset+i][offset+j] = rate_matrices2->getValue()[i].getRate(i, j, 0.0, 1.0);
            }
        }
    }
    
//    // populate rate matrix elements for rm2
//    // off-diagonal blocks of diagonal matrices)
//    for (size_t i = 0; i < num_states_1; i++) {
//        size_t offset = i * num_states_2;
//        for (size_t j = 0; j < num_states_1; j++) {
//            for (size_t k = 0; k < num_states_1; k++) {
//                rate_matrix[offset+i][offset+j] = rate_matrices2->getValue()[i].getRate(i, j, 0.0, 1.0);
//            }
//        }
//    }
    
    // populate rate matrix elements for rm2

    
    // flatten the 2-d rate matrix into a vector
    std::vector<double> all_rates_flat = std::vector<double>( num_states * (num_states - 1), 0.0 );
    size_t k = 0;
    for (size_t i = 0; i < num_states; i++)
    {
        for (size_t j = 0; j < num_states; j++)
        {
            if (i != j)
            {
                all_rates_flat[k++] = rate_matrix[i][j];
            }
        }
    }
    
    // finally set the rates in the actual matrix
    static_cast< RateMatrix_FreeK* >(value)->setTransitionRates( all_rates_flat );
    
    // set the emitted letters for each covartion state in the actual matrix
    std::vector<int> emit(num_states);
    for (int i=0; i<num_states; ++i)
    {
        emit[i] = i % num_states_1;
    }
    static_cast< RateMatrix_FreeK* >(value)->set_emitted_letters( emit);
    
    value->update();
}

void CompositeRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == rate_matrices1)
    {
        rate_matrices1 = static_cast<const TypedDagNode< RbVector<RateGenerator> >* >( newP );
    }
    
    if (oldP == rate_matrices2)
    {
        rate_matrices2 = static_cast<const TypedDagNode< RbVector<RateGenerator> >* >( newP );
    }
    
//    if (oldP == switch_rates)
//    {
//        switch_rates = static_cast<const TypedDagNode<RateGenerator>* >( newP );
//    }
//    
//    if (oldP == clock_rates)
//    {
//        clock_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
//    }
    
}

