//
//  RangeEvolutionRateModifierFunction.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "RangeEvolutionRateModifierFunction.h"

#include "RangeEvolutionRateModifier.h"
#include "TypedDagNode.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

RangeEvolutionRateModifierFunction::RangeEvolutionRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<RbVector<RbVector<double> > >* c, size_t nc) : TypedFunction<CharacterHistoryRateModifier>( new RangeEvolutionRateModifier(nc) ),
gainFactor(gf),
lossFactor(lf),
context_matrix(c)
{

    // add the parameters as parents
    addParameter(gainFactor);
    addParameter(lossFactor);
    addParameter(context_matrix);
    
    update();
}

RangeEvolutionRateModifierFunction::RangeEvolutionRateModifierFunction(const RangeEvolutionRateModifierFunction& m) : TypedFunction<CharacterHistoryRateModifier>( m )
{
    gainFactor = m.gainFactor;
    lossFactor = m.lossFactor;
    context_matrix = m.context_matrix;
}


RangeEvolutionRateModifierFunction::~RangeEvolutionRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}





RangeEvolutionRateModifierFunction* RangeEvolutionRateModifierFunction::clone( void ) const
{
    return new RangeEvolutionRateModifierFunction( *this );
}


void RangeEvolutionRateModifierFunction::update( void )
{
    
    double gf = gainFactor->getValue();
    static_cast<RangeEvolutionRateModifier*>(value)->setGainFactor(gf);
    
    double lf = lossFactor->getValue();
    static_cast<RangeEvolutionRateModifier*>(value)->setLossFactor(lf);
    
    RbVector<RbVector<double> > c = context_matrix->getValue();
    static_cast<RangeEvolutionRateModifier*>(value)->setContextMatrix(c);

}



void RangeEvolutionRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == gainFactor)
    {
        gainFactor = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == lossFactor)
    {
        lossFactor = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == context_matrix)
    {
        context_matrix = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
}
