#ifndef CharacterHistoryDiscrete_H
#define CharacterHistoryDiscrete_H

#include "BranchHistoryDiscrete.h"
#include "CharacterHistory.h"
#include "CharacterEventDiscrete.h"
#include "Cloneable.h"

#include <ostream>
#include <set>
#include <vector>

namespace RevBayesCore {
    
    
    class Tree;
    
    class CharacterHistoryDiscrete : public CharacterHistory {
        
    public:
        
        CharacterHistoryDiscrete(Tree *t, size_t nc, size_t ns, bool rb = false);
        virtual ~CharacterHistoryDiscrete(void);
        
        const BranchHistoryDiscrete&            operator[](size_t i) const;
        
        // public methods
        CharacterHistoryDiscrete*               clone(void) const;
        size_t                                  getNumberStates(void) const;
        CharacterEventDiscrete*                 pickRandomEvent(size_t &bi);
        void                                    setTree(Tree *t);
        
        
    private:
        
        size_t                                  n_states;
        
    };
    
}

#endif

