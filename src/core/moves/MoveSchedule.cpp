#include "MoveSchedule.h"
#include "RbIterator.h"

using namespace RevBayesCore;

MoveSchedule::MoveSchedule(RbVector<Move> *m ) : moves( m )
{
    
}


MoveSchedule::~MoveSchedule()
{
    // we own nothing
}



void MoveSchedule::tune( void )
{
 
    for (RbIterator<Move> it=moves->begin(); it!=moves->end(); ++it)
    {
        it->autoTune();
    }
    
}