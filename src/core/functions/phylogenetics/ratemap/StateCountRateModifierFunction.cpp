//
//  StateCountRateModifierFunction.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "StateCountRateModifierFunction.h"

#include "StateCountRateModifier.h"
#include "TypedDagNode.h"
#include "RbVector.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

StateCountRateModifierFunction::StateCountRateModifierFunction(const TypedDagNode<RbVector<double> >* sf, size_t nc) :
    TypedFunction<CharacterHistoryRateModifier>( new StateCountRateModifier(sf->getValue().size(), nc) ),
    stateFactors(sf)
{
    
    // add the parameters as parents
    
    addParameter(stateFactors);
    
    update();
}

StateCountRateModifierFunction::StateCountRateModifierFunction(const StateCountRateModifierFunction& m) : TypedFunction<CharacterHistoryRateModifier>( m )
{

    stateFactors = m.stateFactors;

}


StateCountRateModifierFunction::~StateCountRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}





StateCountRateModifierFunction* StateCountRateModifierFunction::clone( void ) const
{
    return new StateCountRateModifierFunction( *this );
}


void StateCountRateModifierFunction::update( void )
{
    
    RbVector<double> sf = stateFactors->getValue();
    static_cast<StateCountRateModifier*>(value)->setStateFactors(sf);
    
}



void StateCountRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == stateFactors)
    {
        stateFactors = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
}
