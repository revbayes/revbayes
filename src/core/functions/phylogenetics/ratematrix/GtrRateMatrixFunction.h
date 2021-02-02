#ifndef GtrRateMatrixFunction_H
#define GtrRateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief General Time-Reversible rate matrix function.
     *
     * This function creates the GTR rate matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param er The simplex of exchangeabilities
     * @param bf The simplex of stationary base frequencies
     */
    class GtrRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        GtrRateMatrixFunction(const TypedDagNode< Simplex > *er, const TypedDagNode< Simplex > *bf);
        virtual                                            ~GtrRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        GtrRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                      exchangeability_rates;  //!< A<->C, A<->G, A<->T, C<->G, C<->T, G<->T exchangeabilities
        const TypedDagNode< Simplex >*                      base_frequencies;
        
    };
    
}

#endif
