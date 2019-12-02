#ifndef MinFunction_H
#define MinFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Declaration of the deterministic variable for minimum.
     * The function returns the minimum value within a vector or matrix of doubles
     */
    class MinFunction : public TypedFunction<double> {
        
    public:
        MinFunction(const TypedDagNode< RbVector<double> > * v);
        MinFunction(const TypedDagNode< MatrixReal > * v);
        virtual                                            ~MinFunction(void);                                                          //!< Virtual destructor
        
        // public member functions
        MinFunction*                                        clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        bool                                                matrix;
        const DagNode*                                      vals;
        
    };
    
}


#endif
