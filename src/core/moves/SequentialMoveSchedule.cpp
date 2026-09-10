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
    currentMove( 0 ),
    usedPropOfCurrentMove( 0 )
{
    
    movesPerIteration = 0.0;
    for (RbIterator<Move> it = moves->begin(); it != moves->end(); ++it)
    {
        movesPerIteration += it->getUpdateWeight();
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
    return movesPerIteration;
}


Move& SequentialMoveSchedule::nextMove( std::uint64_t gen )
{
    
    bool found = false;
    do {
        
        double remainingWeight = (*moves)[currentMove].getUpdateWeight() - usedPropOfCurrentMove;
        if ( remainingWeight >= 1.0 )
        {
            usedPropOfCurrentMove++;
            found = true;
        }
        else if ( remainingWeight > 0.0 )
        {
            usedPropOfCurrentMove++;
            RandomNumberGenerator* rng = GLOBAL_RNG;
            found = remainingWeight > rng->uniform01();
        } 
        else
        {
            usedPropOfCurrentMove = 0.0;
            do
            {
                currentMove++;
                if ( currentMove >= moves->size())
                {
                    currentMove = 0;
                }
            } while ( !(*moves)[currentMove].isActive( gen ) );
        }
        
    } while ( !found );
    
    return (*moves)[currentMove];
}
