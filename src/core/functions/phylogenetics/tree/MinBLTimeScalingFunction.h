/**
 * @file
 * This file contains the declaration of the MinBLTimeScaling function, which time-scales
 * a tree using the minimum branch length approach based on a vector of tip ages.
 *
 * @brief Declaration of MinBLTimeScalingFunction
 *
 * @author David Cerny, Laura Mulvey
 * @license GPL version 3
 * @version 1.2.6
 * @since 2024-12-18, version 1.2.6
 *
 */

#ifndef MinBLTimeScalingFunction_H
#define MinBLTimeScalingFunction_H

#include "TypedFunction.h"
#include "Taxon.h"
#include "Tree.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class MinBLTimeScalingFunction : public TypedFunction<Tree> {
        
    public:
        MinBLTimeScalingFunction(const TypedDagNode<Tree> *tr, const TypedDagNode< RbVector<Taxon> > *tx, const TypedDagNode<double> *mbl); //!< Constructor
        virtual                                          ~MinBLTimeScalingFunction(void);     //!< Destructor
        
        // Basic utility functions
        MinBLTimeScalingFunction*                         clone(void) const;                  //!< Clone object
        void                                              update(void);                       //!< Update the function
        
    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP);                 //!< Implementation of swapping parameters
        
    private:
        // members
        const TypedDagNode< Tree >* treeToTimeScale;
        const TypedDagNode< RbVector<Taxon> >* taxonVector;
        const TypedDagNode< double >* minimumBranchLength;
        
    };
    
}

#endif
