#ifndef JcRateMatrixFunction_H
#define JcRateMatrixFunction_H

#include "RateMatrix_JC.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief JC rate matrix function.
     *
     * This function creates the JC rates matrix object by setting the number of states in the constructor.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     * In this case they are actually constant.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class JcRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        JcRateMatrixFunction(size_t ns);
        virtual                                            ~JcRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        JcRateMatrixFunction*                               clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
                
        
    };
    
}

#endif
