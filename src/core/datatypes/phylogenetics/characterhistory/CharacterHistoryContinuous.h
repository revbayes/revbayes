#ifndef CharacterHistoryContinuous_H
#define CharacterHistoryContinuous_H

#include <cstddef>

#include "BranchHistoryContinuous.h"
#include "CharacterEventContinuous.h"
#include "CharacterHistory.h"

namespace RevBayesCore {
    
    
    class Tree;
    
    class CharacterHistoryContinuous : public CharacterHistory {
        
    public:

        CharacterHistoryContinuous(Tree *t, size_t nc, bool rb = false);
        virtual ~CharacterHistoryContinuous(void);
        
        BranchHistoryContinuous&                operator[](size_t i);
        const BranchHistoryContinuous&          operator[](size_t i) const;
        
        CharacterHistoryContinuous*             clone(void) const;
        CharacterEventContinuous*               pickRandomEvent(size_t &bi);
        void                                    setTree(Tree *t);
        
        
    private:
        
    };
    
}

#endif
