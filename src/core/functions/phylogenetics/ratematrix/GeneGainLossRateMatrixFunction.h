/**
 * @file
 * This file contains the declaration of the GeneGainLossRateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the GeneGainLoss rate matrix.
 *
 * @brief Declaration of the GeneGainLossRateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */


#ifndef GeneGainLossRateMatrixFunction_H
#define GeneGainLossRateMatrixFunction_H

#include "AbstractRateMatrix.h"
#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class GeneGainLossRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        GeneGainLossRateMatrixFunction(size_t n, const TypedDagNode<double> *b, const TypedDagNode<double> *d, AbstractRateMatrix::METHOD m);
                
        // public member functions
        GeneGainLossRateMatrixFunction*             clone(void) const;                                                              //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members

        const TypedDagNode<double>*                 birth;
        const TypedDagNode<double>*                 death;

    };
    
}

#endif
