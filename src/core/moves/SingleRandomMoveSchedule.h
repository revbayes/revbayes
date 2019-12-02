#ifndef SingleRandomMoveSchedule_H
#define SingleRandomMoveSchedule_H


#include "Move.h"
#include "MoveSchedule.h"
#include "RbVector.h"

#include <vector>

namespace RevBayesCore {
    
    class SingleRandomMoveSchedule : public MoveSchedule  {
        
    public:
        SingleRandomMoveSchedule(RbVector<Move> *m);                                                                                                                                         //!< Default constructor
        virtual                                        ~SingleRandomMoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        SingleRandomMoveSchedule*                       clone(void) const;
        double                                          getNumberMovesPerIteration(void) const;
        void                                            setNumberMovesPerIteration(double);
        Move&                                           nextMove(unsigned long g);
        
    private:
        
        // Hidden member variables
        double                                          sumOfWeights;
        std::vector<double>                             weights;
    };
    
}

#endif 
