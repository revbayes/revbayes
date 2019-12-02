#include "Kimura81RateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

Kimura81RateMatrixFunction::Kimura81RateMatrixFunction(const TypedDagNode<double> *k1, const TypedDagNode<double> *k2, const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_Kimura81(bf->getValue().size()) ),
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

Kimura81RateMatrixFunction::~Kimura81RateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



Kimura81RateMatrixFunction* Kimura81RateMatrixFunction::clone( void ) const
{
    return new Kimura81RateMatrixFunction( *this );
}


void Kimura81RateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    double k1 = kappa_1->getValue();
    double k2 = kappa_2->getValue();
    const std::vector<double>& f = base_frequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_Kimura81* >(value)->setStationaryFrequencies( f );
    static_cast< RateMatrix_Kimura81* >(value)->setKappa(k1, k2);
    
    value->update();
    
}



void Kimura81RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
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


