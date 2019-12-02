#ifndef LogisticFunction_H
#define LogisticFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Logistic function value of a real number.
     *
     * Compute the logistic function of a real number x (tanh(x) = exp(x) / (1 + exp(x))).
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Nicolas Lartillot)
     * \since Version 1.0, 2015-01-31
     *
     */
    class LogisticFunction : public ContinuousFunction {
        
    public:
        LogisticFunction(const TypedDagNode<double> *a);
        
        LogisticFunction*                   clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         x;
    };
}

#endif
