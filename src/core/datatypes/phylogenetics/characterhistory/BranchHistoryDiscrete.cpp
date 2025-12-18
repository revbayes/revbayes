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

size_t BranchHistoryDiscrete::getMaxObservedState(void) const
{
    size_t max_obs_state = 0;
    
    for (size_t i=0; i<parent_characters.size(); ++i)
    {
        size_t this_state = static_cast<CharacterEventDiscrete*>(parent_characters[i])->getState();
        if ( max_obs_state < this_state )
        {
            max_obs_state = this_state;
        }
    }
    
    for (size_t i=0; i<child_characters.size(); ++i)
    {
        size_t this_state = static_cast<CharacterEventDiscrete*>(child_characters[i])->getState();
        if ( max_obs_state < this_state )
        {
            max_obs_state = this_state;
        }
    }
    
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it = history.begin();
    for (; it!=history.end(); ++it)
    {
        size_t this_state = static_cast<CharacterEventDiscrete*>(*it)->getState();
        if ( max_obs_state < this_state )
        {
            max_obs_state = this_state;
        }
    }
    
    return n_states;
}


size_t BranchHistoryDiscrete::getNumberOfStates(void) const
{
    return n_states;
}


void BranchHistoryDiscrete::setNumberOfStates(size_t n)
{
    n_states = n;
}


CharacterEventDiscrete* BranchHistoryDiscrete::getEvent(size_t i)
{

    return static_cast<CharacterEventDiscrete*>( BranchHistory::getEvent(i) );
}

