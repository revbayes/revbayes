//
#include <iomanip>
#include <cmath>

#include "CharacterEventDiscrete.h"
#include "StateCountRateModifier.h"
#include "RbConstants.h"

using namespace RevBayesCore;

StateCountRateModifier::StateCountRateModifier(size_t ns, size_t nc) : CharacterHistoryRateModifier(ns, nc),
state_factors(  RbVector<double>( 1, num_states) )
{
    ;
}

StateCountRateModifier::StateCountRateModifier(const StateCountRateModifier& g) : CharacterHistoryRateModifier(g)
{
    
    if (&g != this)
    {
        state_factors       = g.state_factors;
    }
}

StateCountRateModifier& StateCountRateModifier::assign(const Assignable &m)
{
    
    const StateCountRateModifier *crm = dynamic_cast<const StateCountRateModifier*>(&m);
    if ( crm != NULL )
    {
        return operator=(*crm);
    }
    else
    {
        throw RbException("Could not assign character history rate modifier.");
    }
}

double StateCountRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age)
{
    
    size_t to_site = newState->getSiteIndex();
    size_t to_state = newState->getState();
    
    double r = 1.0;
    
    // nothing to do if state factor == 1
    if (state_factors[ to_state ] == 1.0) {
        return 1.0;
    }
    
    double n_on = sites_with_states[ to_state ].size();
    if (n_on > 1) {
        double s = state_factors[ to_state ];
        r = std::pow( s, (n_on-1) );
    }
    
    return r;
}

double StateCountRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}

double StateCountRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}


double StateCountRateModifier::computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age)
{
    return 1.0;
}

double StateCountRateModifier::computeSiteRateMultiplier(const TopologyNode& node, unsigned from, unsigned to, unsigned charIdx, double age)
{
    return 1.0;
}


StateCountRateModifier* StateCountRateModifier::clone(void) const
{
    return new StateCountRateModifier(*this);
}

void StateCountRateModifier::update(void)
{
    ; // do nothing
}

void StateCountRateModifier::setStateFactors(RbVector<double> f)
{
    state_factors = f;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const StateCountRateModifier& x)
{
    o << "StateCountRateModifier";
    return o;
}
