#include "DECRateMatrixFunction.h"

#include <cmath>
#include <vector>

#include "RateMatrix_DECRateMatrix.h"
#include "RbMathCombinatorialFunctions.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

DECRateMatrixFunction::DECRateMatrixFunction(   const TypedDagNode< RbVector<RbVector<double> > > *dr,
                                                const TypedDagNode< RbVector<RbVector<double> > > *er,
                                                const TypedDagNode< Simplex > *rs,
                                                bool cs,
                                                bool ex,
//                                                bool os,
                                                bool uc,
                                                size_t mrs) : TypedFunction<RateGenerator>( new RateMatrix_DECRateMatrix((size_t)computeNumStates(er->getValue().size(), mrs, true), er->getValue().size(), cs, ex, uc, mrs) ),
    dispersalRates( dr ),
    extirpationRates( er ),
    rangeSize( rs )
{
    
    // add the rate and frequency parameters as parents
    addParameter( dispersalRates );
    addParameter( extirpationRates );
    addParameter( rangeSize );
    
    update();
}


DECRateMatrixFunction::~DECRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DECRateMatrixFunction* DECRateMatrixFunction::clone( void ) const {
    return new DECRateMatrixFunction( *this );
}

size_t DECRateMatrixFunction::computeNumStates(size_t numAreas, size_t maxRangeSize, bool orderedStates)
{
    if (!orderedStates || maxRangeSize < 1 || maxRangeSize > numAreas)
    {
        return (size_t)pow(2.0, numAreas);
    }
    size_t numStates = 1;
    for (size_t i = 1; i <= maxRangeSize; i++)
    {
        numStates += RbMath::choose(numAreas, i);
    }
    
    return numStates;
}

void DECRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    const RbVector<RbVector<double> >& dr       = dispersalRates->getValue();
    const RbVector<RbVector<double> >& er       = extirpationRates->getValue();
    const std::vector<double>& rs               = rangeSize->getValue();
        
    // set the base frequencies
    static_cast< RateMatrix_DECRateMatrix* >(value)->setDispersalRates(dr);
    static_cast< RateMatrix_DECRateMatrix* >(value)->setExtirpationRates(er);
    static_cast< RateMatrix_DECRateMatrix* >(value)->setRangeSize(rs);
    
    value->update();
}



void DECRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == dispersalRates) {
        dispersalRates = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == extirpationRates) {
        extirpationRates = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == rangeSize) {
        rangeSize = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
}
