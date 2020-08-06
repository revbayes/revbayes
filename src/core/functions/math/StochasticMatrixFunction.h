#ifndef StochasticMatrixFunction_h
#define StochasticMatrixFunction_h

#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"
#include "Simplex.h"

#include <vector>

namespace RevBayesCore {
    
	/**
	 * \brief Create a stochastic matrix.
	 *
	 * This function builds a n x c matrix real with rows that sum to one.
	 *
	 *@param args A vector of simplices. Each simplex is a row in the stochastic matrix.
	 *
	 */
    class StochasticMatrixFunction : public TypedFunction< MatrixReal > {
        
    public:
        StochasticMatrixFunction(const TypedDagNode<RbVector<Simplex > >* &args);
        virtual                                             ~StochasticMatrixFunction(void); //!< Virtual destructor
        
        // public member functions
        StochasticMatrixFunction*                           clone(void) const; //!< Create an independent clone
        void                                                update(void); //! < Update the stochastic matrix when the elements change
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swapping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<Simplex > >*            matrixParams; //! < The vector of simplices that make up the stochasic matrix
        
    };
    
}

#endif /* StochasticMatrixFunction_h */
