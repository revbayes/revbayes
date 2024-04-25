#ifndef SingleRandomMoveSchedule_H
#define SingleRandomMoveSchedule_H


#include <vector>

#include "MoveSchedule.h"

namespace RevBayesCore {
class Move;
template <class valueType> class RbVector;
    
    class SingleRandomMoveSchedule : public MoveSchedule  {
        
    public:
        SingleRandomMoveSchedule(RbVector<Move> *m);                                                                                                                                         //!< Default constructor
        virtual                                        ~SingleRandomMoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        SingleRandomMoveSchedule*                       clone(void) const;
        double                                          getNumberMovesPerIteration(void) const;
        void                                            setNumberMovesPerIteration(double);
        Move&                                           nextMove(std::uint64_t g);
        
    private:
        
        // Hidden member variables
        double                                          sum_of_weights;
        std::vector<double>                             weights;
    };
    
}

#endif 
