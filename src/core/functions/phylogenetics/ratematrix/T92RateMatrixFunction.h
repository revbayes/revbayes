#ifndef T92RateMatrixFunction_H
#define T92RateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief T92 rate matrix function.
     *
     * This function creates the T92 (Tamura 1992) rate matrix object by setting the compound frequency of G+C and the
     * transition-transversion ratio. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param eqGc The compound equilibrium frequency of bases G and C
     * @param tstv The transition-transversion ratio (kappa)
     */
    
    class T92RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        T92RateMatrixFunction(const TypedDagNode<double > *eqGc, const TypedDagNode< double > *tstv);
        virtual                                            ~T92RateMatrixFunction(void);                                                            //!< Virtual destructor
        
        // public member functions
        T92RateMatrixFunction*                              clone(void) const;                                                                      //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< double >*                       equilibriumGc; //!< The compound equilibrium frequency of bases G and C
        const TypedDagNode< double >*                       transitionTransversionRate;
        
    };
    
}

#endif
