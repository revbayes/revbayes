#ifndef RoundFunction_H
#define RoundFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * \brief Round a real number to the next integer.
     *
     * Round a real number to the closest integer.
     * That is, all value x.5 to x.9999 will be rounded up.
     * Otherwise we round down.
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    class RoundFunction : public TypedFunction<long> {
        
    public:
        RoundFunction(const TypedDagNode<double> *a);
        
        RoundFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
