#ifndef MoveSchedule_H
#define MoveSchedule_H

#include "Move.h"
#include "RbVector.h"

#include <vector>

namespace RevBayesCore {
    
    class MoveSchedule {
        
    public:
        virtual                                                ~MoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        virtual MoveSchedule*                                   clone(void) const = 0;
        virtual double                                          getNumberMovesPerIteration(void) const = 0;
        virtual Move&                                           nextMove(unsigned long g) = 0;
        
        // public methods
        void                                                    tune(void);                                                                                     //!< The the moves to achieve better performance.
        
    protected:
        MoveSchedule(RbVector<Move> *m);                                                                                                                  //!< Default constructor
                
        // Hidden member variables
        RbVector<Move>*                                         moves;
        
    };
    
}

#endif
