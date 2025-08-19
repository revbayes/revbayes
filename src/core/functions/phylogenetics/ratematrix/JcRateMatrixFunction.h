#ifndef JcRateMatrixFunction_H
#define JcRateMatrixFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
    
    
    /**
     * @brief JC rate matrix function.
     *
     * This function creates the JC rate matrix object by setting the number of states in the constructor.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     * In this case they are actually constant.
     *
     * @param ns The number of states.
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
