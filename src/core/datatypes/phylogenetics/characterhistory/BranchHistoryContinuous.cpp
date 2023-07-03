#include "BranchHistoryContinuous.h"

#include <vector>

#include "CharacterEventContinuous.h"
#include "Cloneable.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;


BranchHistoryContinuous::BranchHistoryContinuous(size_t nc, size_t idx) : BranchHistory(nc, idx)
{
    
    for (size_t i = 0; i < n_characters; i++)
    {
        parent_characters[i] = new CharacterEventContinuous(i,0.0,0.0);
        child_characters[i] = new CharacterEventContinuous(i,0.0,1.0);
    }
}

BranchHistoryContinuous::~BranchHistoryContinuous(void)
{
    
}


BranchHistoryContinuous* BranchHistoryContinuous::clone(void) const
{
    return new BranchHistoryContinuous(*this);
}


CharacterEventContinuous* BranchHistoryContinuous::getEvent(size_t i)
{
    
    return static_cast<CharacterEventContinuous*>( BranchHistory::getEvent(i) );
}



