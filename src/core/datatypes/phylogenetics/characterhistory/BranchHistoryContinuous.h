#ifndef BranchHistoryContinuous_H
#define BranchHistoryContinuous_H

#include <stddef.h>
#include <set>

#include "BranchHistory.h"
#include "CharacterEventContinuous.h"

namespace RevBayesCore {
    
    class BranchHistoryContinuous : public BranchHistory {
        
    public:

        BranchHistoryContinuous(size_t nc, size_t idx);
        ~BranchHistoryContinuous(void);
        
        // overloaded operators
        bool operator<(const BranchHistoryContinuous&) const;
        
        
        // public methods
        BranchHistoryContinuous*                                        clone(void) const;
        
        // getters
        CharacterEventContinuous*                                       getEvent(size_t i);
        
    protected:
        
        
    };
    
}

#endif

