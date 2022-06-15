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

#ifndef revPoMoThree4RateMatrixFunction_H
#define revPoMoThree4RateMatrixFunction_H

#include "RateMatrix_revPoMoThree4.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "Simplex.h"

#include <vector>

namespace RevBayesCore {

    class revPoMoThree4RateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        revPoMoThree4RateMatrixFunction( const TypedDagNode< Simplex > *bf, const TypedDagNode< RbVector<double> > *ex, const TypedDagNode< RbVector<double> > *fc   );

        virtual                                             ~revPoMoThree4RateMatrixFunction(void);                                                    //!< Virtual destructor

        // public member functions
        revPoMoThree4RateMatrixFunction*                   clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);

    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters

    private:

        // members
        const TypedDagNode< Simplex >*                       pi;
        const TypedDagNode< RbVector<double> >*              rho;
        const TypedDagNode< RbVector<double> >*              phi;


    };

}

#endif
