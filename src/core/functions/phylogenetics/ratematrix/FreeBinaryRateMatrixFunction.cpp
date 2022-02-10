#include "FreeBinaryRateMatrixFunction.h"

#include <vector>

#include "Cloneable.h"
#include "RateMatrix_FreeBinary.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

FreeBinaryRateMatrixFunction::FreeBinaryRateMatrixFunction(const TypedDagNode<RbVector<double> > *tr, bool rescale) : TypedFunction<RateGenerator>( new RateMatrix_FreeBinary(rescale) ),
    transition_rates( tr )
{
    // add the lambda parameter as a parent
    addParameter( transition_rates );
    
    update();
}


FreeBinaryRateMatrixFunction::~FreeBinaryRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



FreeBinaryRateMatrixFunction* FreeBinaryRateMatrixFunction::clone( void ) const
{
    return new FreeBinaryRateMatrixFunction( *this );
}


void FreeBinaryRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    const std::vector<double>& r = transition_rates->getValue();
    
    // set the base frequencies
    static_cast< RateMatrix_FreeBinary* >(value)->setTransitionRates(r);
    
    value->update();
}



void FreeBinaryRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == transition_rates)
    {
        transition_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}

