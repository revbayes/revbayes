#include "CodonSynonymousNonsynonymousHKYRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


CodonSynonymousNonsynonymousHKYRateMatrixFunction::CodonSynonymousNonsynonymousHKYRateMatrixFunction(const TypedDagNode<double> *o, const TypedDagNode<double> *k, const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_CodonSynonymousNonsynonymousHKY() ),
base_frequencies( bf ),
kappa( k ),
omega( o )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    addParameter( kappa );
    addParameter( omega );
    
    update();
}


CodonSynonymousNonsynonymousHKYRateMatrixFunction::~CodonSynonymousNonsynonymousHKYRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CodonSynonymousNonsynonymousHKYRateMatrixFunction* CodonSynonymousNonsynonymousHKYRateMatrixFunction::clone( void ) const
{
    
    return new CodonSynonymousNonsynonymousHKYRateMatrixFunction( *this );
}


void CodonSynonymousNonsynonymousHKYRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    double o = omega->getValue();
    double k = kappa->getValue();
    const std::vector<double>& f = base_frequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_CodonSynonymousNonsynonymousHKY* >(value)->setNucleotideFrequencies( f );
    static_cast< RateMatrix_CodonSynonymousNonsynonymousHKY* >(value)->setKappa( k );
    static_cast< RateMatrix_CodonSynonymousNonsynonymousHKY* >(value)->setOmega( o );
    
    value->update();
}



void CodonSynonymousNonsynonymousHKYRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == omega)
    {
        omega = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}




