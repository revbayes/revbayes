#ifndef BranchHistoryContinuous_H
#define BranchHistoryContinuous_H

#include <cstddef>
#include <set>

#include "BranchHistory.h"
#include "CharacterEventContinuous.h"

namespace RevBayesCore {
    
    class BranchHistoryContinuous : public BranchHistory {
        
    public:

        BranchHistoryContinuous(size_t nc, size_t idx);
        BranchHistoryContinuous(size_t nc, size_t idx, std::set<int> sc);
        BranchHistoryContinuous(const BranchHistoryContinuous& m);
        ~BranchHistoryContinuous(void);
        
        // overloaded operators
        BranchHistoryContinuous& operator=(const BranchHistoryContinuous& bh);
        bool operator<(const BranchHistoryContinuous&) const;
        
        
        // public methods
        BranchHistoryContinuous*                                        clone(void) const;
        
        // getters
        CharacterEventContinuous*                                       getEvent(size_t i);
        
    protected:
        
        
    };
    
}

#endif

