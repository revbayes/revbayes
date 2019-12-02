#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "CharacterEventCompare.h"

#include <iostream>
#include <iomanip>

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

BranchHistoryDiscrete::BranchHistoryDiscrete(size_t nc, size_t ns, size_t idx, std::set<int> sc) : BranchHistory(nc,idx,sc),
    n_states(ns)
{
    
    for (size_t i = 0; i < n_characters; i++)
    {
        parent_characters[i] = new CharacterEventDiscrete(i,0,0.0);
        child_characters[i] = new CharacterEventDiscrete(i,0,1.0);
    }
    
}



BranchHistoryDiscrete::BranchHistoryDiscrete(const BranchHistoryDiscrete& bh) : BranchHistory(bh),
    n_states( bh.n_states )
{
    
}

BranchHistoryDiscrete::~BranchHistoryDiscrete(void)
{
    
}

BranchHistoryDiscrete& BranchHistoryDiscrete::operator=(const BranchHistoryDiscrete &bh)
{
    // delegate to parent class
    BranchHistory::operator=( bh );
    
    if (this != &bh)
    {
        
        n_states                = bh.n_states;

    }
    
    return *this;
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

