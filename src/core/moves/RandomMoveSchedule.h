#ifndef RandomMoveSchedule_H
#define RandomMoveSchedule_H

#include <vector>

#include "MoveSchedule.h"

namespace RevBayesCore {
class Move;
template <class valueType> class RbVector;
    
    class RandomMoveSchedule : public MoveSchedule  {
        
    public:
        RandomMoveSchedule(RbVector<Move> *m);                                                                                                                                         //!< Default constructor
        virtual                                        ~RandomMoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        RandomMoveSchedule*                             clone(void) const;
        double                                          getNumberMovesPerIteration(void) const;
        Move&                                           nextMove(std::uint64_t g);

    private:
        
        // Hidden member variables
        double                                          moves_per_iteration;
        std::vector<double>                             weights;
    };
    
}

#endif
