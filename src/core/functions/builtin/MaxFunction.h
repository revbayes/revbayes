#ifndef MaxFunction_H
#define MaxFunction_H

#include "MatrixReal.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Declaration of the deterministic variable for maximum.
     * The function returns the maximum value within a vector or matrix of doubles
     */
    
    class MaxFunction : public TypedFunction<double> {
        
    public:
        MaxFunction(const TypedDagNode< RbVector<double> > * v);
        MaxFunction(const TypedDagNode< MatrixReal > * v);
        virtual                                            ~MaxFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        MaxFunction*                                        clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swapping parameters
        
    private:
        
        // members
        //!<
        bool                                                matrix;
        const DagNode*                                      vals;
        
    };
    
}


#endif
