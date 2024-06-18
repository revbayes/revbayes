#ifndef CharacterHistory_H
#define CharacterHistory_H

#include <cstddef>
#include <vector>

#include "Cloneable.h"

namespace RevBayesCore {
    
    
    class Tree;
    class BranchHistory;
    class CharacterEvent;
    
    class CharacterHistory : public Cloneable {
        
    public:
        
        virtual ~CharacterHistory(void);
        
        // overloaded operators
        bool                                    operator==(const CharacterHistory &t) const { return false; }
        bool                                    operator!=(const CharacterHistory &t) const { return false; }
        bool                                    operator<(const CharacterHistory &t) const { return false; }
        bool                                    operator<=(const CharacterHistory &t) const { return false; }

        virtual BranchHistory&                  operator[](size_t i);
        virtual const BranchHistory&            operator[](size_t i) const;

        // pure virtual methods
        virtual CharacterHistory*               clone(void) const = 0;
        virtual void                            setTree(const Tree *t) = 0;

        void                                    addEvent(CharacterEvent *e, size_t bi);
        void                                    clear(void);
        const BranchHistory&                    getHistory(size_t n) const;
        size_t                                  getNumberBranches(void) const;
        size_t                                  getNumberEvents(void) const;
        const Tree&                             getTree(void) const;
        bool                                    hasTree(void) const;
        bool                                    hasRootBranch(void) const;
        virtual CharacterEvent*                 pickRandomEvent(size_t &bi);
        void                                    removeEvent(CharacterEvent *e, size_t bi);
        void                                    setHistory(BranchHistory* h, size_t i);
        void                                    setHistories(const std::vector<BranchHistory*>& h);

        
    protected:
        CharacterHistory(const Tree *t, size_t nc, bool rb = false);
        CharacterHistory(const CharacterHistory &ch);
        CharacterHistory&                       operator=(const CharacterHistory &ch);

        const Tree*                             tree;
        std::vector<BranchHistory*>             histories;
        size_t                                  n_branches;
        size_t                                  n_character;
        size_t                                  n_events;
        bool                                    use_root_branch;
        
        
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const CharacterHistory& x);  //!< Overloaded output operator

}

#endif


