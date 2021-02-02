#ifndef CodonSynonymousNonsynonymousHKYRateMatrixFunction_H
#define CodonSynonymousNonsynonymousHKYRateMatrixFunction_H

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
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0.7, 2017-10-16
     *
     */
    class CodonSynonymousNonsynonymousHKYRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        CodonSynonymousNonsynonymousHKYRateMatrixFunction(const TypedDagNode<double> *o, const TypedDagNode<double> *k, const TypedDagNode< Simplex > *bf);
        virtual                                            ~CodonSynonymousNonsynonymousHKYRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        CodonSynonymousNonsynonymousHKYRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone
        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                          base_frequencies;
        const TypedDagNode<double>*                             kappa;
        const TypedDagNode<double>*                             omega;
        
    };
    
}

#endif


