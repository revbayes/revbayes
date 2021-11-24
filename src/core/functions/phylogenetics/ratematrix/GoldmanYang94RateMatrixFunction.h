#ifndef GoldmanYang94RateMatrixFunction_H
#define GoldmanYang94RateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Goldman-Yang (1994) model function.
     *
     * @copyright Copyright 2021-
     * @author Benjamin D Redelings
     * @since 2021-11-24, version 1.0
     */
    class GoldmanYang94RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        GoldmanYang94RateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode<double> *o, const TypedDagNode< Simplex > *pi);
        
        // public member functions
        GoldmanYang94RateMatrixFunction*                    clone(void) const;
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        // members
        const TypedDagNode<double>*                         kappa;
        const TypedDagNode<double>*                         omega;
        const TypedDagNode< Simplex >*                      codon_frequencies;
        
    };
    
}

#endif

