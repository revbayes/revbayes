#ifndef PosteriorPredictiveProbabilityFunction_H
#define PosteriorPredictiveProbabilityFunction_H

#include "RbVector.h"
#include "Simplex.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief Arithmetic PosteriorPredictiveProbability of a vector of real numbers.
     *
     * This function computes the PosteriorPredictiveProbability of a vector of real numbers computing the kth quantile for the vector.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class PosteriorPredictiveProbabilityFunction : public TypedFunction< Simplex > {
        
    public:
        PosteriorPredictiveProbabilityFunction(const TypedDagNode< RbVector<double> > * v, const TypedDagNode<double>* k);
        virtual                                            ~PosteriorPredictiveProbabilityFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        PosteriorPredictiveProbabilityFunction*             clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        const TypedDagNode<double>*                         x;
        
    };
    
}


#endif /* PosteriorPredictiveProbabilityFunction_hpp */

