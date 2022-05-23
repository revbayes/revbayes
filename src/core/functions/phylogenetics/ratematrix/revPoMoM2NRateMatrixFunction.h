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

#ifndef revPoMoM2NRateMatrixFunction_H
#define revPoMoM2NRateMatrixFunction_H

#include "RateMatrix_revPoMoTwo4N.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "Simplex.h"

#include <vector>

namespace RevBayesCore {

    class revPoMoM2NRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        revPoMoM2NRateMatrixFunction( const TypedDagNode< long > *n, const TypedDagNode< RbVector<double> > *m ,  const TypedDagNode< long > *v, const TypedDagNode<double> *g );

        virtual                                             ~revPoMoM2NRateMatrixFunction(void);                                                    //!< Virtual destructor

        // public member functions
        revPoMoM2NRateMatrixFunction*                       clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);

    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters

    private:

        // members
        const TypedDagNode< long >*                          N;
        const TypedDagNode< RbVector<double> >*              mu;
        const TypedDagNode< long >*                          M;
        const TypedDagNode< double >*                        gen;

    };

}

#endif
