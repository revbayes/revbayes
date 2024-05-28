//
#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

#include "CharacterEventDiscrete.h"
#include "RangeEvolutionRateModifier.h"
#include "CharacterHistoryRateModifier.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;

RangeEvolutionRateModifier::RangeEvolutionRateModifier(size_t nc) : CharacterHistoryRateModifier(2, nc),
    gain_factor(0.0),
    loss_factor(0.0),
    context_matrix( std::vector<std::vector<adjacency> >() ),
    forbid_extinction(true),
    is_null_range_absorbing(true)
{
    ;
}

RangeEvolutionRateModifier::RangeEvolutionRateModifier(const RangeEvolutionRateModifier& g) : CharacterHistoryRateModifier(g)
{
    
    if (&g != this)
    {
        gain_factor                = g.gain_factor;
        loss_factor                = g.loss_factor;
        context_matrix             = g.context_matrix;
        forbid_extinction          = g.forbid_extinction;
        is_null_range_absorbing    = g.is_null_range_absorbing;
    }
}

double RangeEvolutionRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age)
{
    
    size_t to_site = newState->getSiteIndex();
    size_t to_state = newState->getState();
    
    double r = 0.0;
    // loss event
    if (to_state == 0)
    {
        if (forbid_extinction && sites_with_states[1].size() == 1)
        {
            // cannot enter the null range (conditions on survival)
            r = 0.0;
            return r;
        }
        else
        {
            r = 1.0;
            return r;
        }
    }
    // gain event
    else if (to_state == 1)
    {
        if (is_null_range_absorbing && sites_with_states[1].size() == 0)
        {
            // cannot leave the null range
            r = 0.0;
            return r;
        }
        else
        {
            for (std::set<size_t>::iterator it = sites_with_states[1].begin(); it != sites_with_states[1].end(); it++)
            {
                size_t from_site = *it;
                r += context_matrix[from_site][to_site].weight;
            }
        }
    }

    return r;
}

double RangeEvolutionRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}

double RangeEvolutionRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}


double RangeEvolutionRateModifier::computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age)
{
    return 1.0;
}

double RangeEvolutionRateModifier::computeSiteRateMultiplier(const TopologyNode& node, unsigned from, unsigned to, unsigned charIdx, double age)
{
    return 1.0;
}


RangeEvolutionRateModifier* RangeEvolutionRateModifier::clone(void) const
{
    return new RangeEvolutionRateModifier(*this);
}

void RangeEvolutionRateModifier::update(void)
{
    ; // do nothing
}

void RangeEvolutionRateModifier::setGainFactor(double f)
{
    gain_factor = f;
}

void RangeEvolutionRateModifier::setLossFactor(double f)
{
    loss_factor = f;
}

void RangeEvolutionRateModifier::setContextMatrix(const RbVector<RbVector<double> >& c)
{
    
    context_matrix = std::vector<std::vector<adjacency> >(this->num_characters);
    
    for (size_t i = 0; i < this->num_characters; i++)
    {
        for (size_t j = 0; j < this->num_characters; j++)
        {
            if (c[i][j] != 0.0)
            {
                adjacency v;
                v.from = i;
                v.to = j;
                v.weight = c[i][j];
                context_matrix[i].push_back(v);
            }
        }
    }
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RangeEvolutionRateModifier& x)
{
    o << "RangeEvolutionRateModifier";
    return o;
}
