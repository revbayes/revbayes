#include <cmath>
#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

#include "CharacterEventDiscrete.h"
#include "HostSwitchRateModifier.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "CharacterHistoryRateModifier.h"
#include "Cloneable.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;

HostSwitchRateModifier::HostSwitchRateModifier(size_t ns, size_t nc) : CharacterHistoryRateModifier(3, nc),
    scale( 1.0 ),
    num_branches( 0 )

{
    ;
}

HostSwitchRateModifier::HostSwitchRateModifier(const HostSwitchRateModifier& g) : CharacterHistoryRateModifier(g)
{
    
    if (&g != this)
    {
        tau = g.tau;
        scale = g.scale;
        distance = g.distance;
        num_branches = g.num_branches;
    }
}

double HostSwitchRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age)
{
    
    // which character will change?
    size_t to_index = newState->getSiteIndex();
    
    // what is the state change?
    size_t from_state = static_cast<CharacterEventDiscrete*>(currState[to_index])->getState();
    size_t to_state = newState->getState();
    
    double r = 1.0;

    
    // loss event (independent of other hosts)
    if (from_state > to_state)
    {
        return 1.0;
    }
    else if (from_state < to_state)
    {
        // rate of leaving 0/1-repertoire equals zero
        size_t num_two = sites_with_states[2].size();
        
        if (num_two == 0) {
            return 0.0;
        }
        
        // gain event
        double scaler_value = scale[ to_state - 1 ];
        
        // if the gain event level's scaling factor equals zero, then there's no effect
        if ( scaler_value == 0.0 )
        {
            return 1.0;
        }
        
        // sum of phylo.distance-scaled rates
        double delta = 0.0;
        size_t n_on = 0;
        for (size_t from_index = 0; from_index < this->num_characters; from_index++)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(currState[from_index])->getState();
            if (s != 0) {
                delta += distance[from_index][to_index];
                n_on += 1;
            } else {
                ; // do nothing
            }
        }
        
        double delta_mean = delta / n_on;
        r = std::pow( delta_mean, -scaler_value);
//        r = std::exp( -scaler_value * delta_mean );
    }
    else {
        throw RbException("Self-transitions not allowed");
    }
    return r;
}

double HostSwitchRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}

double HostSwitchRateModifier::computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age)
{
    std::vector<std::set<size_t> > sites_with_states(num_states);
    for (size_t i = 0; i < currState.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ].insert(i);
    }
    
    return computeRateMultiplier(currState, newState, sites_with_states, age);
}


double HostSwitchRateModifier::computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age)
{
    return 1.0;
}

double HostSwitchRateModifier::computeSiteRateMultiplier(const TopologyNode& node, unsigned from, unsigned to, unsigned charIdx, double age)
{
    return 1.0;
}


HostSwitchRateModifier* HostSwitchRateModifier::clone(void) const
{
    return new HostSwitchRateModifier(*this);
}

void HostSwitchRateModifier::update(void)
{
    ; // do nothing
}

void HostSwitchRateModifier::setTree(const RevBayesCore::Tree &t)
{
    tau = t;
    num_branches = tau.getNumberOfNodes() - 1;
    distance = *RevBayesCore::TreeUtilities::getDistanceMatrix ( tau );
    
    double max_distance = 2 * tau.getRoot().getAge();
    double sum_distance = 0.0;
    for (size_t i = 0; i < distance.size(); i++) {
        for (size_t j = i; j < distance[i].size(); j++) {
            distance[i][j] /= max_distance;
            distance[j][i] /= max_distance;
            sum_distance += distance[i][j];
        }
    }
    
    size_t n_tips = distance.size();
    size_t n_pairs = (n_tips*n_tips-n_tips) / 2;
    double mean_distance = sum_distance / n_pairs;
    for (size_t i = 0; i < distance.size(); i++) {
        for (size_t j = i; j < distance[i].size(); j++) {
            distance[i][j] /= mean_distance;
            distance[j][i] /= mean_distance;
        }
    }
    

}

void HostSwitchRateModifier::setScale(const std::vector<double>& s)
{
    scale = s;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const HostSwitchRateModifier& x)
{
    o << "HostSwitchRateModifier";
    return o;
}
