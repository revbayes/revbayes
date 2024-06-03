//
//  StateCountRateModifier.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef StateCountRateModifier_h
#define StateCountRateModifier_h

#include "CharacterHistoryRateModifier.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <set>
#include <string>
#include <vector>


namespace RevBayesCore
{
    
    
    
    class StateCountRateModifier : public CharacterHistoryRateModifier
    {
    public:
        StateCountRateModifier(size_t ns, size_t nc);
        StateCountRateModifier(const StateCountRateModifier& g);
        
        
        
        double                              computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, unsigned currState, unsigned newState, unsigned charIdx=0, double age=0.0);
        
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age=0.0);
        

        void                                setStateFactors(RbVector<double> f);
        
        
        void                                update(void);
        StateCountRateModifier*             clone(void) const;
        
    protected:
        
        
    private:
        
        RbVector<double>                    state_factors;
        
    };
    
    std::ostream& operator<<(std::ostream& o, const StateCountRateModifier& x);
}

#endif /* StateCountRateModifier_h */
