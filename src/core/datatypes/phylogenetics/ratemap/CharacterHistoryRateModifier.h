#ifndef CharacterHistoryRateModifier_H
#define CharacterHistoryRateModifier_H

#include <set>
#include <vector>
#include "Cloneable.h"

namespace RevBayesCore
{
    class CharacterEvent;
    class CharacterEventDiscrete;
    class CharacterHistoryRateModifier : public Cloneable
    {
    public:
        CharacterHistoryRateModifier(size_t ns, size_t nc);
        CharacterHistoryRateModifier(const CharacterHistoryRateModifier& g);
        
        bool                                operator==(const CharacterHistoryRateModifier &rm) const { return this == &rm; }
        bool                                operator!=(const CharacterHistoryRateModifier &rm) const { return !operator==(rm); }
        bool                                operator<(const CharacterHistoryRateModifier &rm) const { return this < &rm; }
        bool                                operator<=(const CharacterHistoryRateModifier &rm) const { return operator<(rm) || operator==(rm); }
        
        virtual double                      computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age=0.0) = 0;
        virtual double                      computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<size_t> counts, double age=0.0);
        virtual double                      computeRateMultiplier(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, std::vector<std::set<size_t> > sites_with_states, double age=0.0);
        virtual std::set<size_t>            getAffectedSites(CharacterEventDiscrete* newState) const;

        virtual void                        update(void) = 0;
        CharacterHistoryRateModifier*       clone( void ) const = 0;

    protected:
        size_t                              num_states;
        size_t                              num_characters;

    private:
        std::set<size_t>                    all_sites;

    };
    std::ostream& operator<<(std::ostream& o, const CharacterHistoryRateModifier& x);                                         //!< Overloaded output operator
    std::ostream& operator<<(std::ostream& o, const std::vector<CharacterHistoryRateModifier*>& x);                                         //!< Overloaded output operator
}

#endif
