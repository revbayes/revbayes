#ifndef OrderedRateMatrixFunction_H
#define OrderedRateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief Ordered rate matrix function.
     *
     * This function creates the ordered rate matrix object by setting the number of states, the rate of gains, and
     * the rate of losses. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param n The number of states
     * @param l The rate of gains
     * @param m The number of losses
     * @param allow_zero_state Should state '0' be allowed? (May not be appropriate for some counts)
     */
    
    class OrderedRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        OrderedRateMatrixFunction(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *l, const TypedDagNode<double> *m, bool allow_zero_state, bool rescale, std::string method);
        
        virtual                                     ~OrderedRateMatrixFunction(void);                                               //!< Virtual destructor
        
        // public member functions
        OrderedRateMatrixFunction*                  clone(void) const;                                                              //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode<double>*                 lambda;
        const TypedDagNode<double>*                 mu;
        bool                                        zero;

    };
    
}

#endif
