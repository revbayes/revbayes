#ifndef FreeSymmetricRateMatrixFunction_H
#define FreeSymmetricRateMatrixFunction_H

#include <iosfwd>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Free symmetric rate matrix function.
     *
     * This function creates the free (symmetric) rate matrix object by setting the substitution rates.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param tr The vector of substitution rates, whose number of elements equals half the number of the off-diagonal entries of the rate matrix.
     * @param r Should the rates be rescaled so that the average rate equals 1?
     * @param method Matrix exponentiation method.
     */
    class FreeSymmetricRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        FreeSymmetricRateMatrixFunction(const TypedDagNode< RbVector<double> > *tr , bool r, std::string method);
        virtual                                            ~FreeSymmetricRateMatrixFunction(void);                                      //!< Virtual destructor
        
        // public member functions
        FreeSymmetricRateMatrixFunction*                    clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             transition_rates;
        
    };
    
}

#endif
