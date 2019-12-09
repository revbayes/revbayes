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

#ifndef PoMoTwoRateMatrixFunction_H
#define PoMoTwoRateMatrixFunction_H

#include "RateMatrix_PoMoTwo.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "Simplex.h"

#include <vector>

namespace RevBayesCore {

    class PoMoTwoRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        PoMoTwoRateMatrixFunction(const TypedDagNode< long > *ps, const TypedDagNode< RbVector<double> > *rho, const TypedDagNode< Simplex > *pi );

        virtual                                             ~PoMoTwoRateMatrixFunction(void);                                                    //!< Virtual destructor

        // public member functions
        PoMoTwoRateMatrixFunction*                          clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);

    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters

    private:

        // members
        const TypedDagNode< long >*                          populationSize;
        const TypedDagNode< RbVector<double> >*              exchangeabilities;
        const TypedDagNode< Simplex >*                       equilibriumFrequencies;

    };

}

#endif
