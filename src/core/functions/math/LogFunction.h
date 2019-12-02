#ifndef LogFunction_H
#define LogFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * \brief Logarithm of a to base b.
     *
     * The logarithm of value 'a' taken to base 'b'.
     * Useful application are for base 2, e and 10.
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    class LogFunction : public ContinuousFunction {
        
    public:
        LogFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b);
        
        LogFunction*                        clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         a;
        const TypedDagNode<double>*         base;
    };
}

#endif
