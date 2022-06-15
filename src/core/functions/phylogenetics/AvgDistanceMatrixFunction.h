/**
 * @file
 * This file contains the declaration of the AvgDistanceMatrix function, which computes
 * the average distance matrix from a vector of distance matrices.
 *
 * @brief Declaration of the AvgDistanceMatrixFunction
 *
 * @author David Cerny
 * @license GPL version 3
 * @version 1.1
 * @since 2021-10-19, version 1.1
 *
 */

#ifndef AvgDistanceMatrixFunction_H
#define AvgDistanceMatrixFunction_H

#include "TypedFunction.h"
#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class AvgDistanceMatrixFunction : public TypedFunction< AverageDistanceMatrix > {
        
    public:
        AvgDistanceMatrixFunction(const TypedDagNode< RbVector<DistanceMatrix> >* matvect);   //!< Constructor for the unweighted version
        AvgDistanceMatrixFunction(const TypedDagNode< RbVector<DistanceMatrix> >* matvect, const TypedDagNode< RbVector<double> >* weights);                                                                                      //!< Constructor for the weighted version
        virtual                                          ~AvgDistanceMatrixFunction(void);    //!< Destructor
        
        // Basic utility functions
        AvgDistanceMatrixFunction*                        clone(void) const;                  //!< Clone object
        void                                              update(void);                       //!< Clone the function
        
    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP);                 //!< Implementation of swapping parameters
        
    private:
        // members
        const TypedDagNode< RbVector<DistanceMatrix> >* matrixVector;
        const TypedDagNode< RbVector<double> >* weightVector;
        
    };
    
}

#endif
