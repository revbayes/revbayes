#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbIterator.h"
#include "SequentialMoveSchedule.h"
#include "Move.h"
#include "MoveSchedule.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

SequentialMoveSchedule::SequentialMoveSchedule(RbVector<Move> *s) : MoveSchedule( s ),
    current_move( 0 ),
    used_prop_of_current_move( 0 )
{
    
    moves_per_iteration = 0.0;
    for (RbIterator<Move> it = moves->begin(); it != moves->end(); ++it)
    {
        moves_per_iteration += it->getUpdateWeight();
    }
}


SequentialMoveSchedule::~SequentialMoveSchedule()
{
    // we own nothing
}


SequentialMoveSchedule* SequentialMoveSchedule::clone( void ) const
{
    return new SequentialMoveSchedule(*this);
}

double SequentialMoveSchedule::getNumberMovesPerIteration( void ) const
{
    return moves_per_iteration;
}


Move& SequentialMoveSchedule::nextMove( std::uint64_t gen )
{
    
    bool found = false;
    do {
        
        double remainingWeight = (*moves)[current_move].getUpdateWeight() - used_prop_of_current_move;
        if ( remainingWeight >= 1.0 )
        {
            used_prop_of_current_move++;
            found = true;
        }
        else if ( remainingWeight > 0.0 )
        {
            used_prop_of_current_move++;
            RandomNumberGenerator* rng = GLOBAL_RNG;
            found = remainingWeight > rng->uniform01();
        } 
        else
        {
            used_prop_of_current_move = 0.0;
            do
            {
                current_move++;
                if ( current_move >= moves->size())
                {
                    current_move = 0;
                }
            } while ( !(*moves)[current_move].isActive( gen ) );
        }
        
    } while ( !found );
    
    return (*moves)[current_move];
}
