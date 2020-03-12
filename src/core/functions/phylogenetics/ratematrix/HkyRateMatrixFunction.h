#ifndef HkyRateMatrixFunction_H
#define HkyRateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief HKY rate matrix function.
     *
     * This function creates the HKY rates matrix object by setting the transition-transversion parameter kappa
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param k The transition-transversion ratio (kappa)
     * @param bf The simplex of base frequencies
     *
     */
    class HkyRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        HkyRateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode< Simplex > *bf);
        virtual                                            ~HkyRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        HkyRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                      base_frequencies;
        const TypedDagNode<double>*                         kappa;
        
    };
    
}

#endif
