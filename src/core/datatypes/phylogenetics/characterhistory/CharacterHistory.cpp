#include "CharacterHistory.h"

#include <cassert>
#include <math.h>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "BranchHistory.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;


CharacterHistory::CharacterHistory(Tree *t, size_t nc, bool rb ) :
    tree( t ),
    histories(),
    n_branches(),
    n_character( nc ),
    n_events( 0 ),
    use_root_branch( rb )
{
    
}


CharacterHistory::CharacterHistory(const CharacterHistory &ch) :
    tree( ch.tree ),
    histories(),
    n_branches( ch.n_branches ),
    n_character( ch.n_character ),
    n_events( ch.n_events ),
    use_root_branch( ch.use_root_branch )
{
    
    for (size_t i=0; i<ch.histories.size(); ++i)
    {
        BranchHistory *bh = ch.histories[i];
        histories.push_back( bh->clone() );
    }
    
}


/**
 * Dummy destructor.
 */
CharacterHistory::~CharacterHistory()
{
 
    for (size_t i=0; i<histories.size(); ++i)
    {
        BranchHistory *bh = histories[i];
        delete bh;
    }
    
}


CharacterHistory& CharacterHistory::operator=(const CharacterHistory &ch)
{
    
    if ( this != &ch )
    {
    
        // first, free the memory
        for (size_t i=0; i<histories.size(); ++i)
        {
            BranchHistory *bh = histories[i];
            delete bh;
        }
        
        // now, assign the values
        tree            = ch.tree;
        n_branches      = ch.n_branches;
        n_character     = ch.n_character;
        n_events        = ch.n_events;
        use_root_branch = ch.use_root_branch;

        histories.clear();
        for (size_t i=0; i<ch.histories.size(); ++i)
        {
            BranchHistory *bh = ch.histories[i];
            histories.push_back( bh->clone() );
        }

    }
    
    return *this;
}


/**
 * Index operator.
 */
const BranchHistory& CharacterHistory::operator[](size_t i) const
{
    if ( i > histories.size() )
    {
        throw RbException("Index out of bounds in character history.");
    }
    
    return *(histories[i]);
}


/**
 * Add an event onto the tree.
 * Here we will also
 */
void CharacterHistory::addEvent( CharacterEvent *e, size_t branch_index)
{
    
    assert( branch_index < histories.size() );
    
    BranchHistory *bh = histories[ branch_index ];
    bh->addEvent( e );
    
    ++n_events;
    
    size_t counted_events = 0;
    for ( size_t i=0; i<histories.size(); ++i )
    {
        counted_events += histories[ i ]->getNumberEvents();
    }
    assert( counted_events == n_events );
}


/**
 * Get the number of branches of the tree.
 */
size_t CharacterHistory::getNumberBranches( void ) const
{
    return n_branches;
}



/**
 * Get the number of total events over the tree.
 */
size_t CharacterHistory::getNumberEvents( void ) const
{
    return n_events;
}


bool CharacterHistory::hasRootBranch( void ) const
{
    return use_root_branch;
}

/**
 * Pick a random event.
 */
CharacterEvent* CharacterHistory::pickRandomEvent( size_t &branch_index )
{
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    // randomly pick an index for the event
    size_t event_index = size_t( std::floor(n_events * rng->uniform01()) );
    for (size_t i=0; i<n_branches; ++i)
    {
        BranchHistory* bh = histories[i];
        if ( bh->getNumberEvents() > event_index )
        {
            branch_index = i;
            return bh->getEvent(event_index);
        }
        else
        {
            event_index -= bh->getNumberEvents();
        }
        
    }
    
    throw RbException("Could not pick a random event because I think there are more events than there actually are.");
}



/**
 * Remove the given event from our histories.
 */
void CharacterHistory::removeEvent( CharacterEvent *e, size_t branch_index)
{
    size_t counted_events_before = 0;
    for ( size_t i=0; i<histories.size(); ++i )
    {
        counted_events_before += histories[ i ]->getNumberEvents();
    }
    assert( counted_events_before == n_events );
    
    BranchHistory *bh = histories[ branch_index ];
    bh->removeEvent( e );
    
    --n_events;
    
    size_t counted_events = 0;
    for ( size_t i=0; i<histories.size(); ++i )
    {
        counted_events += histories[ i ]->getNumberEvents();
    }
    assert( counted_events == n_events );
}
