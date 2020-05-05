/**
 * @file
 * This file contains the declaration of the OrderedRateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the Ordered rate matrix.
 *
 * @brief Declaration of the OrderedRateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */


#ifndef OrderedRateMatrixFunction_H
#define OrderedRateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class OrderedRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        OrderedRateMatrixFunction(const TypedDagNode<long> *n, const TypedDagNode<double> *l, const TypedDagNode<double> *m, bool allow_zero_state);
        
        virtual                                     ~OrderedRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        OrderedRateMatrixFunction*                  clone(void) const;                                                              //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode<long>*                   n;
        const TypedDagNode<double>*                 lambda;
        const TypedDagNode<double>*                 mu;

    };
    
}

#endif
