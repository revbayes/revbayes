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
    
    // Loss event
    if (from_state > to_state)
    {
        // Braga et al. (2020) allow only repertoires holding at least one actual
        // host (2), so any transition to such a state should have rate zero.
        if (from_state == 2 && sites_with_states[2].size() == 1)
        {
            return 0.0;
        }
        else 
        {
            // The loss rate is otherwise unaffected by the rest of the repertoire
            return 1.0;
        }
    }
    // Gain event
    else if (from_state < to_state)
    {
        // A repertoire without any actual hosts (2s) is unreachable from any 
        // valid repertoire. This ensures that invalid repertoires sampled at the root
        // propagate will be invalid at the leaves, and are rejected by the MCMC
        size_t num_two = sites_with_states[2].size();
        if (num_two == 0) {
            return 0.0;
        }
        
        // Read the current value of beta
        double beta = scale[ to_state - 1 ];
        
        // If the gain event level's scaling factor equals zero, then there's no effect
        if ( beta == 0.0 )
        {
            return 1.0;
        }
        
        // The phylogenetic factor enters as the average normalized phylogenetic
        // distance to the the set of potentials and/or actual hosts.
        // Normalization of the host phylogeny is performed by setTree().
        double delta = 0.0;
        size_t n_on = 0;
        for (size_t from_index = 0; from_index < this->num_characters; from_index++)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(currState[from_index])->getState();
            // For a 0 -> 1 transition all potential (1) and actual (2) hosts affect the rate
            // For a 1 -> 2 transition only the actual (2) hosts affect the rate
            bool include = (to_state == 2) ? (s == 2) : (s != 0);
            if (include) {
                delta += distance[from_index][to_index];
                n_on += 1;
            }
        }

        double delta_mean = delta / n_on;
        r = std::exp( -beta * delta_mean );
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
