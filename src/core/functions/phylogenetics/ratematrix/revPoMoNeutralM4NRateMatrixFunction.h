/*
 * @file
 * This file contains the declaration of the reversible Pomo2 rate matrix function class.
 * This class is derived from the function class and is used to compute the rate matrix of a general time reversible Markov chain.
 *
 * @brief Declaration of the Pomo2 rate matrix function.
 *
 * @date Last modified: 19/07/2019
 * @author Rui Borges
 * @license GPL version 3
 * @version 1.0
 * @interface Function
 *
 * $Id$
 */

#ifndef revPoMoNeutralM4NRateMatrixFunction_H
#define revPoMoNeutralM4NRateMatrixFunction_H

#include "RateMatrix_revPoMoNeutralM4N.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "Simplex.h"

#include <vector>
#include <cstdint>

namespace RevBayesCore {

    class revPoMoNeutralM4NRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        revPoMoNeutralM4NRateMatrixFunction( const TypedDagNode< std::int64_t > *n, const TypedDagNode< std::int64_t > *m,  const TypedDagNode< Simplex > *bf, const TypedDagNode< RbVector<double> > *ex  );

        virtual                                             ~revPoMoNeutralM4NRateMatrixFunction(void);                                                    //!< Virtual destructor

        // public member functions
        revPoMoNeutralM4NRateMatrixFunction*                clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);

    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters

    private:

        std::int64_t                                                 computeNumStates( std::int64_t mi );

        // members
        const TypedDagNode< std::int64_t >*                          N;
        const TypedDagNode< std::int64_t >*                          M;
        const TypedDagNode< Simplex >*                       pi;
        const TypedDagNode< RbVector<double> >*              rho;

    };

}

#endif
