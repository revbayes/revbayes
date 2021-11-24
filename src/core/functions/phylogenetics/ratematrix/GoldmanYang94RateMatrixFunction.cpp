#include "GoldmanYang94RateMatrixFunction.h"

#include <vector>

#include "Cloneable.h"
#include "RateMatrix_GoldmanYang94.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


GoldmanYang94RateMatrixFunction::GoldmanYang94RateMatrixFunction(const TypedDagNode<double> *k,
                                                                 const TypedDagNode<double> *o,
                                                                 const TypedDagNode< Simplex > *pi)
    : TypedFunction<RateGenerator>( new RateMatrix_GoldmanYang94() ),
      kappa( k ),
      omega( o ),
      codon_frequencies( pi )
{
    addParameter( kappa );
    addParameter( omega );
    addParameter( codon_frequencies );
    
    update();
}


GoldmanYang94RateMatrixFunction* GoldmanYang94RateMatrixFunction::clone( void ) const
{
    return new GoldmanYang94RateMatrixFunction( *this );
}


void GoldmanYang94RateMatrixFunction::update( void )
{
    // get the argument values for the function
    double k = kappa->getValue();
    double o = omega->getValue();
    const std::vector<double>& f = codon_frequencies->getValue();
    
    // set the argument values for the function
    static_cast< RateMatrix_GoldmanYang94* >(value)->setCodonFrequencies( f );
    static_cast< RateMatrix_GoldmanYang94* >(value)->setOmega( o );
    static_cast< RateMatrix_GoldmanYang94* >(value)->setKappa( k );
    
    value->update();
}



void GoldmanYang94RateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == codon_frequencies)
    {
        codon_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    if (oldP == omega)
    {
        omega = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}



