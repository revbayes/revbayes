//
//  SiteRateModifier.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/23/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef SiteRateModifier_h
#define SiteRateModifier_h

#include "CharacterHistoryRateModifier.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <set>
#include <cstdint>
#include <string>
#include <vector>


namespace RevBayesCore
{
    
    
    
    class SiteRateModifier : public CharacterHistoryRateModifier
    {
    public:
        SiteRateModifier(size_t ns, size_t nc);
        SiteRateModifier(const SiteRateModifier& g);
        
        
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, CharacterEvent* currState, CharacterEvent* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, unsigned currState, unsigned newState, unsigned charIdx=0, double age=0.0);
        
        void                                setRateMultipliers(const RbVector<RbVector<double> >& rm);
        void                                setEventClasses(const RbVector<RbVector<std::int64_t> >& ec);
        void                                setSiteClasses(const RbVector<std::int64_t>& sc);
        
        void                                update(void);
        SiteRateModifier*                   clone(void) const;
        
    protected:
        
        
    private:
        RbVector<RbVector<double> >         rate_multipliers;
        RbVector<RbVector<std::int64_t> >            event_classes;
        RbVector<std::int64_t>                       site_classes;
        
        size_t                              num_event_classes;
        size_t                              num_site_classes;
    };
    
    std::ostream& operator<<(std::ostream& o, const SiteRateModifier& x);
}

#endif /* SiteRateModifier_h */
