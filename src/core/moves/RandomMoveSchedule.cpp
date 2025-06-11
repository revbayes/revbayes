#include "RandomMoveSchedule.h"

#include <cstddef>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbIterator.h"
#include "Move.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

RandomMoveSchedule::RandomMoveSchedule(RbVector<Move> *s) : MoveSchedule( s )
{
    
    moves_per_iteration = 0.0;
    for (RbIterator<Move> it = moves->begin(); it != moves->end(); ++it)
    {
        moves_per_iteration += it->getUpdateWeight();
        weights.push_back( it->getUpdateWeight() );
    }
}


RandomMoveSchedule::~RandomMoveSchedule()
{
    // we own nothing
}


RandomMoveSchedule* RandomMoveSchedule::clone( void ) const
{
    return new RandomMoveSchedule(*this);
}


double RandomMoveSchedule::getNumberMovesPerIteration( void ) const
{
    return moves_per_iteration;
}


Move& RandomMoveSchedule::nextMove( std::uint64_t gen )
{
    
    moves_per_iteration = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if ( (*moves)[i].isActive( gen ) )
        {
            moves_per_iteration += weights[i];
        }
    }
    
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    double u = moves_per_iteration * rng->uniform01();
    
    size_t index = 0;
    // only if the move is inactive or the weight of the move is smaller than u
    while ( !(*moves)[index].isActive(gen) || weights[index] <= u )
    {
        // check if this move is active
        // if not, then we just subtract the weight of this move
        if ( (*moves)[index].isActive( gen ) )
        {
            u -= weights[index];
        }
        ++index;
    }

    if (index >= moves->size())
    {
        index = moves->size() - 1;
    }
        
    return (*moves)[index];
}
