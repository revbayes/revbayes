#include <cstddef>
#include <vector>

#include "BranchHistoryDiscrete.h"
#include "CharacterHistoryDiscrete.h"
#include "RbException.h"
#include "Tree.h"
#include "CharacterEventDiscrete.h"
#include "CharacterHistory.h"
#include "Cloneable.h"

namespace RevBayesCore { class BranchHistory; }

using namespace RevBayesCore;


CharacterHistoryDiscrete::CharacterHistoryDiscrete(Tree *t, size_t nc, size_t ns, bool rb ) : CharacterHistory(t, nc, rb),
    n_states( ns )
{
    
    if ( tree != NULL )
    {
        n_branches = tree->getNumberOfNodes() - 1 + ( use_root_branch ? 1 : 0 );
        // create a branch history object for each branch
        for (size_t i=0; i<n_branches; ++i)
        {
            histories.push_back( new BranchHistoryDiscrete(n_character, n_states, i) );
        }
        
    }
    
}



/**
 * Dummy destructor.
 */
CharacterHistoryDiscrete::~CharacterHistoryDiscrete()
{
    
}


/**
 * Index operator.
 */
const BranchHistoryDiscrete& CharacterHistoryDiscrete::operator[](size_t i) const
{
    if ( i > histories.size() )
    {
        throw RbException("Index out of bounds in character history.");
    }
    
    return static_cast<const BranchHistoryDiscrete &>(*histories[i]);
}


/**
 * Create a deep copy.
 */
CharacterHistoryDiscrete* CharacterHistoryDiscrete::clone( void ) const
{
    
    return new CharacterHistoryDiscrete( *this );
}


/**
 * Get the number of states for this character.
 */
size_t CharacterHistoryDiscrete::getNumberStates( void ) const
{
    return n_states;
}


/**
 * Pick a random event.
 */
CharacterEventDiscrete* CharacterHistoryDiscrete::pickRandomEvent( size_t &branch_index )
{
    
    return static_cast<CharacterEventDiscrete *>( CharacterHistory::pickRandomEvent(branch_index) );
}


/**
 * Set the new tree.
 */
void CharacterHistoryDiscrete::setTree(Tree *t)
{
    
    tree = t;
    if ( tree != NULL )
    {
        histories.clear();
        n_events = 0;
        n_branches = tree->getNumberOfNodes() - 1  + (use_root_branch ? 1 : 0);
        
        // create a branch history object for each branch
        for (size_t i=0; i<n_branches; ++i)
        {
            histories.push_back( new BranchHistoryDiscrete(n_character, n_states, i) );
        }
        
    }
    
}






