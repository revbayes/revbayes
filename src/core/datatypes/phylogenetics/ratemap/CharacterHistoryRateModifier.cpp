//
//  CharacterHistoryRateModifier.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 8/8/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include <iostream>
#include <cstddef>
#include <set>
#include <vector>

#include "CharacterHistoryRateModifier.h"
#include "RbException.h"
#include "Assignable.h"

namespace RevBayesCore { class CharacterEvent; }
namespace RevBayesCore { class CharacterEventDiscrete; }

using namespace RevBayesCore;

CharacterHistoryRateModifier::CharacterHistoryRateModifier(size_t ns, size_t nc) :
    num_states(ns),
    num_characters(nc)
{
    for (size_t i = 0; i < num_characters; i++) {
        all_sites.insert(i);
    }
}


CharacterHistoryRateModifier::CharacterHistoryRateModifier(const CharacterHistoryRateModifier& g)
{
    
    if (&g != this)
    {
        all_sites = g.all_sites;
        num_states = g.num_states;
        num_characters = g.num_characters;
    }
}

CharacterHistoryRateModifier& CharacterHistoryRateModifier::assign(const Assignable &m)
{
    const CharacterHistoryRateModifier *crm = dynamic_cast<const CharacterHistoryRateModifier*>(&m);
    if ( crm != NULL )
    {
        return operator=(*crm);
    }
    else
    {
        throw RbException("Could not assign character history rate modifier.");
    }
}


double CharacterHistoryRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age)
{
    return computeRateMultiplier(currState, newState);
}

double CharacterHistoryRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age)
{
    
    return computeRateMultiplier(currState, newState);
}


std::set<size_t> CharacterHistoryRateModifier::getAffectedSites(CharacterEventDiscrete* newState) const
{
    return all_sites;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<CharacterHistoryRateModifier*>& x)
{
    o << " ";
    
    return o;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const CharacterHistoryRateModifier& x)
{
    o << " ";
    
    return o;
}
