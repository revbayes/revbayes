//
//  HostSwitchRateModifier.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/24/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef HostSwitchRateModifier_h
#define HostSwitchRateModifier_h


#include "CharacterHistoryRateModifier.h"
#include "DistanceMatrix.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

#include <set>
#include <string>
#include <vector>


namespace RevBayesCore
{
    
    
    
    class HostSwitchRateModifier : public CharacterHistoryRateModifier
    {
    public:
        HostSwitchRateModifier(size_t ns, size_t nc);
        HostSwitchRateModifier(const HostSwitchRateModifier& g);
        
        
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, unsigned currState, unsigned newState, unsigned charIdx=0, double age=0.0);
        
        void                                setTree(const Tree& t);
        void                                setScale(const std::vector<double>& s);
        
        void                                update(void);
        HostSwitchRateModifier*   clone(void) const;
        
    protected:
        
        
    private:
        DistanceMatrix                      distance;
        Tree                                tau;
        std::vector<double>                 scale;
        size_t                              num_branches;
        std::vector<std::vector<double> >   node_phylo_adjacencies;
        std::vector<std::vector<double> >   node_phylo_distances;
    };
    
    std::ostream& operator<<(std::ostream& o, const HostSwitchRateModifier& x);
}

#endif /* HostSwitchRateModifier_h */
