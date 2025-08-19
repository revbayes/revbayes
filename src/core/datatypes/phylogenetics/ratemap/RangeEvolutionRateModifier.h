//
//  RangeEvolutionRateModifier.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef RangeEvolutionRateModifier_h
#define RangeEvolutionRateModifier_h

#include "CharacterHistoryRateModifier.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <set>
#include <string>
#include <vector>


namespace RevBayesCore
{
    
    
    
    class RangeEvolutionRateModifier : public CharacterHistoryRateModifier
    {
    public:
        RangeEvolutionRateModifier(size_t nc);
        RangeEvolutionRateModifier(const RangeEvolutionRateModifier& g);
        
        

        double                              computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, unsigned currState, unsigned newState, unsigned charIdx=0, double age=0.0);
        
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age=0.0);
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age=0.0);
        
        
        void                                setContextMatrix(const RbVector<RbVector<double> >& c);
        void                                setGainFactor(double f);
        void                                setLossFactor(double f);
        
        
        void                                update(void);
        RangeEvolutionRateModifier*         clone(void) const;
        
    protected:
        
        
    private:
        void                                initializeContexts(RbVector<RbVector<double> > c);
        
        struct adjacency {
            size_t                          from;
            size_t                          to;
            double                          weight;
        };
        
        std::vector<std::vector<adjacency> > context_matrix;
        double                               gain_factor;
        double                               loss_factor;
        
        bool                                 forbid_extinction;
        bool                                 is_null_range_absorbing;
    
    };
    
    std::ostream& operator<<(std::ostream& o, const RangeEvolutionRateModifier& x);
}

#endif /* RangeEvolutionRateModifier_h */
