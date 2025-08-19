#ifndef InfiniteSitesRateMatrixFunction_H
#define InfiniteSitesRateMatrixFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
    
    
    /**
     * @brief Infinite dites rate matrix function.
     *
     * This function creates the infinite sites rates matrix object by setting the number of states in the constructor.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     * In this case they are actually constant.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-04-30
     *
     */
    class InfiniteSitesRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        InfiniteSitesRateMatrixFunction(size_t ns);
        virtual                                            ~InfiniteSitesRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        InfiniteSitesRateMatrixFunction*                    clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
        
    };
    
}

#endif
