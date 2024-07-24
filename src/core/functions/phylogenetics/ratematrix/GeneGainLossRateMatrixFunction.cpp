#include "GeneGainLossRateMatrixFunction.h"

#include "RateMatrix_GeneGainLoss.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

GeneGainLossRateMatrixFunction::GeneGainLossRateMatrixFunction(size_t n, const TypedDagNode<double> *b, const TypedDagNode<double> *d, AbstractRateMatrix::METHOD m) : TypedFunction<RateGenerator>( new RateMatrix_GeneGainLoss(n, m) ),
    birth( b ),
    death( d )
{

    addParameter( birth );
    addParameter( death );
    
    update();
}



GeneGainLossRateMatrixFunction* GeneGainLossRateMatrixFunction::clone( void ) const
{
    return new GeneGainLossRateMatrixFunction( *this );
}


void GeneGainLossRateMatrixFunction::update( void )
{
    double b = birth->getValue();
    double d = death->getValue();

    static_cast< RateMatrix_GeneGainLoss* >(value)->setBirth( b );
    static_cast< RateMatrix_GeneGainLoss* >(value)->setDeath( d );

    static_cast< RateMatrix_GeneGainLoss* >(value)->update();
    
}



void GeneGainLossRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == birth)
    {
        birth = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == death)
    {
        death = static_cast<const TypedDagNode<double>* >( newP );
    }
        
}



