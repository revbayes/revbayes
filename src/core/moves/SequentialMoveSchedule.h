/**
 * @file
 * This file contains the declaration of Move, which is the abstract
 * base class for moves used in the mcmc inference.
 *
 * @brief Declaration of the base class Move
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-17, version 1.0
 *
 * $Id$
 */

#ifndef SequentialMoveSchedule_H
#define SequentialMoveSchedule_H

#include "Move.h"
#include "MoveSchedule.h"
#include "RbVector.h"

#include <vector>
#include <cstdint>

namespace RevBayesCore {
    
    class SequentialMoveSchedule : public MoveSchedule  {
        
    public:
        SequentialMoveSchedule(RbVector<Move> *m);                                                                                                                                         //!< Default constructor
        virtual                                        ~SequentialMoveSchedule(void);                                                                             //!< Destructor
        
        // pure virtual public methods
        SequentialMoveSchedule*                         clone(void) const;
        double                                          getNumberMovesPerIteration(void) const;
        Move&                                           nextMove(std::uint64_t g);
        
    private:
        
        // Hidden member variables
        double                                          movesPerIteration;
        size_t                                          currentMove;
        double                                          usedPropOfCurrentMove;
    };
    
}

#endif
