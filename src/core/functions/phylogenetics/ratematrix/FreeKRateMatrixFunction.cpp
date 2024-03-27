#include "FreeKRateMatrixFunction.h"

#include <math.h>
#include <cstddef>
#include <vector>

#include "RateMatrix_FreeK.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

FreeKRateMatrixFunction::FreeKRateMatrixFunction(const TypedDagNode< RbVector<double> > *trf, bool rescale, std::string method) : TypedFunction<RateGenerator>( new RateMatrix_FreeK( 0.5+sqrt(0.25+trf->getValue().size() ), rescale, method ) ),
    transition_rates( NULL ),
    transition_rates_flat( trf )
{
    
    // add the rate and frequency parameters as parents
    addParameter( transition_rates_flat );
    
    update();
}

FreeKRateMatrixFunction::FreeKRateMatrixFunction(const TypedDagNode< RbVector<RbVector<double> > >*tr, bool rescale, std::string method) : TypedFunction<RateGenerator>( new RateMatrix_FreeK( tr->getValue().size(), rescale, method ) ),
    transition_rates( tr ),
    transition_rates_flat( NULL )
{
    
    // add the rate and frequency parameters as parents
    addParameter( transition_rates );
    
    update();
}

FreeKRateMatrixFunction::~FreeKRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



FreeKRateMatrixFunction* FreeKRateMatrixFunction::clone( void ) const
{
    return new FreeKRateMatrixFunction( *this );
}


void FreeKRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    if (transition_rates != NULL && transition_rates_flat == NULL)
    {
        const RbVector<RbVector<double> >& r = transition_rates->getValue();
        
        size_t n = r.size();
        std::vector<double> r_flat( n * (n-1) );
        size_t k = 0;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                if (i != j)
                {
                    r_flat[k++] = r[i][j];
                }
                
            }
            
        }

        // set the flattened rates
        static_cast< RateMatrix_FreeK* >(value)->setTransitionRates(r_flat);
    }
    else if (transition_rates == NULL && transition_rates_flat != NULL)
    {
        const std::vector<double>& r = transition_rates_flat->getValue();
        
        // set the flattened rates
        static_cast< RateMatrix_FreeK* >(value)->setTransitionRates(r);
    }

    value->update();
}



void FreeKRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == transition_rates)
    {
        transition_rates = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == transition_rates_flat)
    {
        transition_rates_flat = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

