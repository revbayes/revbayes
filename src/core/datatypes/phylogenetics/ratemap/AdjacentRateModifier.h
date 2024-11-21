#ifndef AdjacentRateModifier_H
#define AdjacentRateModifier_H

#include "CharacterHistoryRateModifier.h"
#include "MatrixReal.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <set>
#include <cstdint>
#include <string>
#include <vector>


namespace RevBayesCore
{
    
    class AdjacentRateModifier : public CharacterHistoryRateModifier
    {
    public:
        AdjacentRateModifier(size_t ns, size_t nc);
        AdjacentRateModifier(const AdjacentRateModifier& g);
        
        
        double                              computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, CharacterEventDiscrete* currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeSiteRateMultiplier(const TopologyNode& node, unsigned currState, unsigned newState, unsigned charIdx=0, double age=0.0);
        
        double                              computeRateMultiplierUsingWidth(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        double                              computeRateMultiplierUsingMatrix(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0);
        virtual std::set<size_t>            getAffectedSites(CharacterEventDiscrete* newState) const;

        
        void                                setContextMatrix(const RbVector<RbVector<std::int64_t> >& c);
        void                                setContextMatrix(const MatrixReal& c);
        void                                setGainFactor(double f);
        void                                setLossFactor(double f);
        void                                setWidth(size_t w);
        
        
        void                                update(void);
        AdjacentRateModifier*               clone(void) const;
        
    protected:
       
        
    private:
        void                                initializeContexts(RbVector<RbVector<std::int64_t> > c);
        
        struct adjacency {
            size_t                          from;
            size_t                          to;
            double                          weight;
        };
        
        std::vector<std::vector<adjacency> > context_matrix;
        size_t                              width;
        double                              gain_factor;
        double                              loss_factor;
        std::vector<double>                 exp_gain_factors;
        std::vector<double>                 exp_loss_factors;
        std::vector<std::set<size_t> >      context_set;
        
        std::string                         context_type;
    };
    
    std::ostream& operator<<(std::ostream& o, const AdjacentRateModifier& x);
}

#endif /* defined(__rb_mlandis__AdjacentRateModifier__) */
