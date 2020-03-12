#ifndef K80RateMatrixFunction_H
#define K80RateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief K80 rate matrix function.
     *
     * This function creates the Kimura80 rate matrix object by setting the transition-transversion parameter kappa.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param k The transition-transversion ratio (kappa)
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
