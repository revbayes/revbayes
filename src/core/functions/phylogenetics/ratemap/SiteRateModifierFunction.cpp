//
//  SiteRateModifierFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "SiteRateModifierFunction.h"

#include "SiteRateModifier.h"
#include "TypedDagNode.h"
#include "RbVector.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

SiteRateModifierFunction::SiteRateModifierFunction(const TypedDagNode<RbVector<RbVector<double> > >* rm,
                                                   const TypedDagNode<RbVector<RbVector<std::int64_t> > >* ec,
                                                   const TypedDagNode<RbVector<std::int64_t> >* sc) :
    TypedFunction<CharacterHistoryRateModifier>( new SiteRateModifier(ec->getValue().size(), sc->getValue().size() ) ),
    rate_multipliers(rm),
    event_classes(ec),
    site_classes(sc)
{
    // add the parameters as parents
    addParameter(rate_multipliers);
    addParameter(event_classes);
    addParameter(site_classes);
    
    update();
}

SiteRateModifierFunction::SiteRateModifierFunction(const SiteRateModifierFunction& m) : TypedFunction<CharacterHistoryRateModifier>( m )
{
    rate_multipliers = m.rate_multipliers;
    event_classes    = m.event_classes;
    site_classes     = m.site_classes;
}


SiteRateModifierFunction::~SiteRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}





SiteRateModifierFunction* SiteRateModifierFunction::clone( void ) const
{
    return new SiteRateModifierFunction( *this );
}


void SiteRateModifierFunction::update( void )
{
    
    // get values to update
    const RbVector<RbVector<double> >& rm = rate_multipliers->getValue();
    const RbVector<RbVector<std::int64_t> >& ec    = event_classes->getValue();
    const RbVector<std::int64_t>& sc               = site_classes->getValue();
    
    // apply updates to SiteRateModifier
    static_cast<SiteRateModifier*>(value)->setRateMultipliers(rm);
    static_cast<SiteRateModifier*>(value)->setEventClasses(ec);
    static_cast<SiteRateModifier*>(value)->setSiteClasses(sc);
    
}



void SiteRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == rate_multipliers)
    {
        rate_multipliers = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == event_classes)
    {
        event_classes = static_cast<const TypedDagNode<RbVector<RbVector<std::int64_t> > >* >( newP );
    }
    else if (oldP == site_classes)
    {
        site_classes = static_cast<const TypedDagNode<RbVector<std::int64_t> >* >( newP );
    }
}
