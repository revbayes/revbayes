#ifndef RandomMoveSchedule_H
#define RandomMoveSchedule_H

#include "Move.h"
#include "MoveSchedule.h"
#include "RbVector.h"

#include <vector>

namespace RevBayesCore {
    
    class RandomMoveSchedule : public MoveSchedule  {
        
    public:
        RandomMoveSchedule(RbVector<Move> *m);                                                                                                                                         //!< Default constructor
        virtual                                        ~RandomMoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        RandomMoveSchedule*                             clone(void) const;
        double                                          getNumberMovesPerIteration(void) const;
        Move&                                           nextMove(unsigned long g);

    private:
        
        // Hidden member variables
        double                                          movesPerIteration;
        std::vector<double>                             weights;
    };
    
}

#endif
