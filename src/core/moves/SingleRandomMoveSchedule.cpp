#include "SingleRandomMoveSchedule.h"

#include <cstddef>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbIterator.h"
#include "Move.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

SingleRandomMoveSchedule::SingleRandomMoveSchedule(RbVector<Move> *s) : MoveSchedule( s )
{
    
    sumOfWeights = 0.0;
    for (RbIterator<Move> it = moves->begin(); it != moves->end(); ++it)
    {
        sumOfWeights+= it->getUpdateWeight();
        weights.push_back( it->getUpdateWeight() );
    }
}


SingleRandomMoveSchedule::~SingleRandomMoveSchedule()
{
    // we own nothing
}


SingleRandomMoveSchedule* SingleRandomMoveSchedule::clone( void ) const
{

    return new SingleRandomMoveSchedule(*this);
}

double SingleRandomMoveSchedule::getNumberMovesPerIteration( void ) const
{
    return 1.0;
}


Move& SingleRandomMoveSchedule::nextMove( std::uint64_t gen )
{
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    double u = sumOfWeights * rng->uniform01();
    
    size_t index = 0;
    while ( weights[index] < u || !(*moves)[index].isActive( gen ) )
    {
        u -= weights[index];
        ++index;
    }
    
    return (*moves)[index];
}
