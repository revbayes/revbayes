#include <vector>

#include "RateGenerator.h"
#include "TamuraNeiRateMatrixFunction.h"
#include "TypedFunction.h"
#include "RateMatrix_TamuraNei.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Default constructor.
 *
 * This function takes three inputs:
 * @param k1 The ratio of A<->G transitions to transversions
 * @param k2 The ratio of C<->T transitions to transversions
 * @param bf The simplex of base frequencies
 */

TamuraNeiRateMatrixFunction::TamuraNeiRateMatrixFunction(const TypedDagNode<double> *k1, const TypedDagNode<double> *k2, const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_TamuraNei() ),
    kappa_1( k1 ),
    kappa_2( k2 ),
    base_frequencies( bf )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    addParameter( kappa_1 );
    addParameter( kappa_2 );
    
    update();
}

TamuraNeiRateMatrixFunction::~TamuraNeiRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TamuraNeiRateMatrixFunction* TamuraNeiRateMatrixFunction::clone( void ) const
{
    return new TamuraNeiRateMatrixFunction( *this );
}


void TamuraNeiRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    double k1 = kappa_1->getValue();
    double k2 = kappa_2->getValue();
    const std::vector<double>& f = base_frequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_TamuraNei* >(value)->setStationaryFrequencies( f );
    static_cast< RateMatrix_TamuraNei* >(value)->setKappa(k1, k2);
    
    value->update();
    
}



void TamuraNeiRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode<Simplex>* >( newP );
    }
    
    if (oldP == kappa_1)
    {
        kappa_1 = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == kappa_2)
    {
        kappa_2 = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}


