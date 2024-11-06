#include <cstddef>
#include <cstdint>
#include <ostream>
#include <set>
#include <vector>

#include "CharacterEventDiscrete.h"
#include "SiteRateModifier.h"
#include "CharacterHistoryRateModifier.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbIterator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;

SiteRateModifier::SiteRateModifier(size_t ns, size_t nc) : CharacterHistoryRateModifier(ns, nc),
    rate_multipliers( RbVector<RbVector<double> >() ),
    event_classes( RbVector<RbVector<std::int64_t> >() ),
    site_classes( RbVector<std::int64_t>() ),
    num_event_classes( 0 ),
    num_site_classes( 0 )
{
    ;
}

SiteRateModifier::SiteRateModifier(const SiteRateModifier& g) : CharacterHistoryRateModifier(g)
{
    
    if (&g != this)
    {
        rate_multipliers = g.rate_multipliers;
        event_classes = g.event_classes;
        site_classes = g.site_classes;
        num_event_classes = g.num_event_classes;
        num_site_classes = g.num_site_classes;
    }
}

double SiteRateModifier::computeRateMultiplier(std::vector<CharacterEvent *> currState, CharacterEventDiscrete* newState, double age)
{
    
    // the new event acts on which site?
    size_t new_site = newState->getSiteIndex();
    
    // which site_class does that site belong to?
    size_t new_site_class = site_classes[new_site];
    
    // the new event is what event type?
    size_t old_state = static_cast<CharacterEventDiscrete*>(currState[new_site])->getState();
    size_t new_state = newState->getState();
    
    // which event_class does that event belong to?
    size_t new_event_class = event_classes[old_state][new_state];
    
    // pick (site_class,event_class) from rate_multipliers
    double f = rate_multipliers[new_site_class][new_event_class] * num_event_classes;
    
    
    return f;
}

double SiteRateModifier::computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age)
{
    return 1.0;
}

double SiteRateModifier::computeSiteRateMultiplier(const TopologyNode& node, unsigned from, unsigned to, unsigned charIdx, double age)
{
    return 1.0;
}


SiteRateModifier* SiteRateModifier::clone(void) const
{
    return new SiteRateModifier(*this);
}

void SiteRateModifier::update(void)
{
    ; // do nothing
}

void SiteRateModifier::setRateMultipliers(const RbVector<RbVector<double> >& rm)
{
    rate_multipliers = rm;
}

void SiteRateModifier::setEventClasses(const RbVector<RbVector<std::int64_t> >& ec)
{
    event_classes = ec;
    std::set<size_t> s;
    for (size_t i = 0; i < event_classes.size(); i++)
    {
        s.insert( event_classes[i].begin(), event_classes[i].end() );
    }
    num_event_classes = s.size();
}

void SiteRateModifier::setSiteClasses(const RbVector<std::int64_t>& sc)
{
    site_classes = sc;
    std::set<size_t> s( site_classes.begin(), site_classes.end() );
    num_site_classes = s.size();
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const SiteRateModifier& x)
{
    o << "SiteRateModifier";
    return o;
}
