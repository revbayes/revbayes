#ifndef F81RateMatrixFunction_H
#define F81RateMatrixFunction_H

#include "RateMatrix_F81.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief F81 rate matrix function.
     *
     * This function creates the F81 rates matrix object by setting the base frequencies.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class F81RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        F81RateMatrixFunction(const TypedDagNode< Simplex > *bf);
        virtual                                            ~F81RateMatrixFunction(void);                                                            //!< Virtual destructor
        
        // public member functions
        F81RateMatrixFunction*                              clone(void) const;                                                                      //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                      base_frequencies;
        
    };
    
}

#endif
