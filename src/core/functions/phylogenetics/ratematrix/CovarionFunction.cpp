#include "CovarionFunction.h"

#include <cstddef>
#include <vector>

#include "RateMatrix_FreeK.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

CovarionFunction::CovarionFunction(bool r) : TypedFunction<RateGenerator>( NULL ),
    rescale( r ),
    rate_matrices( NULL ),
    rate_scalars( NULL ),
    switch_rates( NULL )
{
    
}


CovarionFunction::~CovarionFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


CovarionFunction* CovarionFunction::clone( void ) const
{
    return new CovarionFunction( *this );
}


void CovarionFunction::setRateMatrices(const TypedDagNode< RbVector<RateGenerator> > *rm)
{
    rate_matrices = rm;
    addParameter( rate_matrices );
    
    if ( rate_scalars != NULL && switch_rates != NULL && rate_matrices != NULL )
    {
        size_t num_states = rate_scalars->getValue().size() * rate_matrices->getValue()[0].getNumberOfStates();
        value = new RateMatrix_FreeK( num_states, rescale );
        update();
    }
}


void CovarionFunction::setRateScalars(const TypedDagNode< RbVector<double> > *rs)
{
    rate_scalars = rs;
    
    addParameter( rate_scalars );
    
    if ( rate_scalars != NULL && switch_rates != NULL && rate_matrices != NULL )
    {
        size_t num_states = rate_scalars->getValue().size() * rate_matrices->getValue()[0].getNumberOfStates();
        value = new RateMatrix_FreeK( num_states, rescale );
        update();
    }
}



void CovarionFunction::setSwitchRates(const TypedDagNode< RbVector<RbVector<double> > > *sr)
{
    switch_rates = sr;
    
    addParameter( switch_rates );
    
    if ( rate_scalars != NULL && switch_rates != NULL && rate_matrices != NULL )
    {
        size_t num_states = rate_scalars->getValue().size() * rate_matrices->getValue()[0].getNumberOfStates();
        value = new RateMatrix_FreeK( num_states, rescale );
        update();
    }
}


void CovarionFunction::update( void )
{
    
    
    
    size_t num_org_states = rate_matrices->getValue()[0].getNumberOfStates();
    size_t num_categories = rate_scalars->getValue().size();
    
    // set up a 2-d matrix to hold the rates for the combined observed and hidden states
    size_t num_states = num_org_states * num_categories;
    std::vector< std::vector<double> > rate_matrix = std::vector< std::vector<double> >( num_states, std::vector<double>( num_states, 0.0 ) );
    for (size_t initial_category = 0; initial_category < num_categories; initial_category++)
    {
        for (size_t initial_org_state = 0; initial_org_state < num_org_states; initial_org_state++)
        {
            size_t i = (initial_category * num_org_states) + initial_org_state;
            
            for (size_t final_category = 0; final_category < num_categories; final_category++)
            {
                for (size_t final_org_state = 0; final_org_state < num_org_states; final_org_state++)
                {
                    size_t j = (final_category * num_org_states) + final_org_state;
                    
                    if (initial_category == final_category)
                    {
                        const RateGenerator &rg = rate_matrices->getValue()[initial_category];
                        rate_matrix[i][j] = rg.getRate(initial_org_state, final_org_state, 0.0, 1.0) * rate_scalars->getValue()[initial_category];
                    }
                    else if (initial_org_state == final_org_state)
                    {
                        rate_matrix[i][j] = switch_rates->getValue()[initial_category][final_category];
                    }
                }
            }
        }
    }
    
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

    // set the emitted letters for each covarion state in the actual matrix
    std::vector<int> emit(num_states);
    for (int i=0; i<num_states; ++i)
    {
        emit[i] = i % num_org_states;
    }
    static_cast< RateMatrix_FreeK* >(value)->set_emitted_letters( emit);

    value->update();
}


void CovarionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == rate_matrices)
    {
        rate_matrices = static_cast<const TypedDagNode< RbVector<RateGenerator> >* >( newP );
    }
    else if (oldP == rate_scalars)
    {
        rate_scalars = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == switch_rates)
    {
        switch_rates = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    
}
