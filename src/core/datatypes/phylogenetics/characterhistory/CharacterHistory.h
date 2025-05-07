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
        
        virtual const BranchHistory&             operator[](size_t i) const;

        // pure virtual methods
        virtual CharacterHistory*               clone(void) const = 0;
        virtual void                            setTree(Tree *t) = 0;

        void                                    addEvent(CharacterEvent *e, size_t bi);
        size_t                                  getNumberBranches(void) const;
        size_t                                  getNumberEvents(void) const;
        bool                                    hasRootBranch(void) const;
        virtual CharacterEvent*                 pickRandomEvent(size_t &bi);
        void                                    removeEvent(CharacterEvent *e, size_t bi);
        
        
    protected:
        CharacterHistory(Tree *t, size_t nc, bool rb = false);
        CharacterHistory(const CharacterHistory &ch);
        CharacterHistory&                       operator=(const CharacterHistory &ch);

        Tree*                                   tree;
        std::vector<BranchHistory*>             histories;
        size_t                                  n_branches;
        size_t                                  n_character;
        size_t                                  n_events;
        bool                                    use_root_branch;
        
        
    };
    
}

#endif


