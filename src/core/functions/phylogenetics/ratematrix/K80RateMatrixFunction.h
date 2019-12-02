#ifndef K80RateMatrixFunction_H
#define K80RateMatrixFunction_H

#include "RateMatrix_Kimura80.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief K80 rate matrix function.
     *
     * This function creates the K80 rates matrix object by setting the transition-transversion parameter kappa.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-11-18
     *
     */
    class K80RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        K80RateMatrixFunction(const TypedDagNode<double> *k);
        virtual                                            ~K80RateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        K80RateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                         kappa;
        
    };
    
}

#endif
