//
//  HostSwitchRateModifierFunction.cpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/24/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "HostSwitchRateModifierFunction.h"

#include <vector>

#include "HostSwitchRateModifier.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

HostSwitchRateModifierFunction::HostSwitchRateModifierFunction(const TypedDagNode<Tree>* t,
                                                               const TypedDagNode<RbVector<double> >* s) :
TypedFunction<CharacterHistoryRateModifier>( new HostSwitchRateModifier(3, t->getValue().getNumberOfTips() ) ),
tau(t),
scale(s)
{
    // add the parameters as parents
    addParameter(tau);
    addParameter(scale);
    
    update();
}

HostSwitchRateModifierFunction::HostSwitchRateModifierFunction(const HostSwitchRateModifierFunction& m) : TypedFunction<CharacterHistoryRateModifier>( m )
{
    tau  = m.tau;
    scale = m.scale;
}


HostSwitchRateModifierFunction::~HostSwitchRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


HostSwitchRateModifierFunction* HostSwitchRateModifierFunction::clone( void ) const
{
    return new HostSwitchRateModifierFunction( *this );
}


void HostSwitchRateModifierFunction::update( void )
{
    
    // get values to update
    const Tree& t = tau->getValue();
    const std::vector<double>& s = scale->getValue();
    
    // apply updates to HostSwitchRateModifier
    static_cast<HostSwitchRateModifier*>(value)->setTree(t);
    static_cast<HostSwitchRateModifier*>(value)->setScale(s);

    
}



void HostSwitchRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == scale)
    {
        scale = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
}
