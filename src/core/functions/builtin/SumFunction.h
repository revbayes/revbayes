#ifndef SumFunction_H
#define SumFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Declaration of the deterministic variable for minimum.
     * The function returns the sum of a vector or matrix of doubles.
     */
    class SumFunction : public TypedFunction<double> {
        
    public:
        SumFunction(const TypedDagNode< RbVector<double> > * v);
        SumFunction(const TypedDagNode< MatrixReal > * v);
        virtual                                            ~SumFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        SumFunction*                                        clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swapping parameters
        
    private:
        
        // members
        bool                                                matrix;
        const DagNode*                                      vals;
        
    };
    
}


#endif
