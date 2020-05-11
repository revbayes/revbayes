#include "WagRateMatrixFunction.h"

#include <vector>

#include "Cloneable.h"
#include "RateMatrix_Wag.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

WagRateMatrixFunction::WagRateMatrixFunction( const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_Wag() ),
base_frequencies( bf )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    
    update();
}

WagRateMatrixFunction::~WagRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



WagRateMatrixFunction* WagRateMatrixFunction::clone( void ) const
{
    return new WagRateMatrixFunction( *this );
}


void WagRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    const std::vector<double>& f = base_frequencies->getValue();
    
    // set the base frequencies
    static_cast< RateMatrix_Wag* >(value)->setStationaryFrequencies( f );
    
    value->update();
    
}



void WagRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
}


