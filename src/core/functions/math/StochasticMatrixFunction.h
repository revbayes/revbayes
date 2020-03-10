#ifndef StochasticMatrixFunction_h
#define StochasticMatrixFunction_h

#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"
#include "Simplex.h"

#include <vector>

namespace RevBayesCore {
    
    class StochasticMatrixFunction : public TypedFunction< MatrixReal > {
        
    public:
        StochasticMatrixFunction(const TypedDagNode<RbVector<Simplex > >* &args);
        virtual                                             ~StochasticMatrixFunction(void);                                                      //!< Virtual destructor
        
        // public member functions
        StochasticMatrixFunction*                           clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<Simplex > >*            matrixParams;
        
    };
    
}

#endif /* StochasticMatrixFunction_h */
