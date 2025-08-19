#ifndef BranchHistory_H
#define BranchHistory_H

#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

#include "CharacterEventCompare.h"
#include "Cloneable.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class CharacterEvent;

    class BranchHistory : public Cloneable {

    public:
        ~BranchHistory(void);
        
        // overloaded operators
        BranchHistory& operator=(const BranchHistory& bh);
        //bool operator()(const BranchHistory&) const;
        bool operator<(const BranchHistory&) const;


        // pure virtual methods
        virtual BranchHistory*                                          clone(void) const = 0;

        // public methods
        bool                                                            areEventTimesValid(const TopologyNode &node) const;
        void                                                            addEvent(CharacterEvent* evt);
        void                                                            clearEvents(void);
        void                                                            clearEvents(const std::set<size_t>& clearSet);

        // getters
        const size_t                                                    getBranchIndex(void) const;
        virtual CharacterEvent*                                         getEvent(size_t i);
        const size_t                                                    getNumberCharacters(void) const;
        const size_t                                                    getNumberEvents(void) const;
        std::vector<CharacterEvent*>&                                   getParentCharacters(void);
        const std::vector<CharacterEvent*>&                             getParentCharacters(void) const;
        std::vector<CharacterEvent*>&                                   getChildCharacters(void);
        const std::vector<CharacterEvent*>&                             getChildCharacters(void) const;
        std::multiset<CharacterEvent*,CharacterEventCompare>&           getHistory(void);
        const std::multiset<CharacterEvent*,CharacterEventCompare>&     getHistory(void) const;

        void                                                            print(const TopologyNode* nd=NULL) const;

        // setters
        void                                                            setParentCharacters(const std::set<CharacterEvent*>& s);
        void                                                            setParentCharacters(const std::vector<CharacterEvent*>& s);
        void                                                            setChildCharacters(const std::set<CharacterEvent*>& s);
        void                                                            setChildCharacters(const std::vector<CharacterEvent*>& s);
        void                                                            setHistory(const std::set<CharacterEvent*,CharacterEventCompare>& s);
        void                                                            setHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& s);

        void                                                            removeEvent(CharacterEvent* evt);
        void                                                            updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<CharacterEvent*>& parentSet, const std::set<CharacterEvent*>& childSet, const std::set<size_t>& indexSet);
        void                                                            updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<size_t>& indexSet);
        void                                                            updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet);

        

    protected:
        //BranchHistory(void);
        BranchHistory(size_t nc, size_t idx);
        BranchHistory(size_t nc, size_t idx, std::set<int> sc);
        BranchHistory(const BranchHistory& m);

        // container/element arguments
        size_t                                                          n_characters;
        size_t                                                          branch_index;

        // containers
        mutable std::multiset<CharacterEvent*,CharacterEventCompare>    history;
        mutable std::vector<CharacterEvent*>                            parent_characters;
        mutable std::vector<CharacterEvent*>                            child_characters;


    };

    std::ostream& operator<<(std::ostream& o, const BranchHistory& x);
}

#endif
