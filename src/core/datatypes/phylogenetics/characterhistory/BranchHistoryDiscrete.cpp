#include "BranchHistoryDiscrete.h"

#include <vector>

#include "CharacterEventDiscrete.h"
#include "Cloneable.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;



BranchHistoryDiscrete::BranchHistoryDiscrete(size_t nc, size_t ns, size_t idx) : BranchHistory(nc,idx),
    n_states(ns)
{
    
    for (size_t i = 0; i < n_characters; i++)
    {
        parent_characters[i] = new CharacterEventDiscrete(i,0,0.0);
        child_characters[i] = new CharacterEventDiscrete(i,0,1.0);
    }
}

BranchHistoryDiscrete::~BranchHistoryDiscrete(void)
{
    
}


BranchHistoryDiscrete* BranchHistoryDiscrete::clone(void) const
{
    return new BranchHistoryDiscrete(*this);
}


const size_t BranchHistoryDiscrete::getNumberStates(void) const
{
    return n_states;
}


CharacterEventDiscrete* BranchHistoryDiscrete::getEvent(size_t i)
{

    return static_cast<CharacterEventDiscrete*>( BranchHistory::getEvent(i) );
}

