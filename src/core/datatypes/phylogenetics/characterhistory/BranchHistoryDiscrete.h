#ifndef BranchHistoryDiscrete_H
#define BranchHistoryDiscrete_H

#include <cstddef>
#include <ostream>
#include <set>

#include "BranchHistory.h"
#include "CharacterEventDiscrete.h"

namespace RevBayesCore {
    
    class BranchHistoryDiscrete : public BranchHistory {
        
    public:
        //BranchHistory(void);
        BranchHistoryDiscrete(size_t nc, size_t ns, size_t idx);
        BranchHistoryDiscrete(size_t nc, size_t ns, size_t idx, std::set<int> sc);
        BranchHistoryDiscrete(const BranchHistoryDiscrete& m);
        ~BranchHistoryDiscrete(void);
        
        // overloaded operators
        BranchHistoryDiscrete& operator=(const BranchHistoryDiscrete& bh);
        
        // public methods
        BranchHistoryDiscrete*                                          clone(void) const;
        
        // getters
        CharacterEventDiscrete*                                         getEvent(size_t i);
        const size_t                                                    getNumberStates(void) const;
        
        
    protected:
        
        // container/element arguments
        size_t                                                          n_states;
        
        
    };
    
    std::ostream& operator<<(std::ostream& o, const BranchHistoryDiscrete& x);
}

#endif

