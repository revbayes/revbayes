#ifndef Kimura81RateMatrixFunction_H
#define Kimura81RateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Kimura81 rate matrix function.
     *
     * This function creates the Kimura81 rate matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param k1 The ratio of transitions (A<->G, C<->T) to A<->C, G<->T transversions
     * @param k2 The ratio of A<->T, C<->G transversions to A<->C, G<->T transversions
     * @param bf The simplex of stationary base frequencies
     */
    class Kimura81RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        Kimura81RateMatrixFunction(const TypedDagNode<double> *k1, const TypedDagNode<double> *k2, const TypedDagNode< Simplex > *bf);
        virtual                                            ~Kimura81RateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        Kimura81RateMatrixFunction*                         clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode<double>*                         kappa_1; //!< The ratio of transitions (A<->G, C<->T) to A<->C, G<->T transversions
        const TypedDagNode<double>*                         kappa_2; //!< The ratio of A<->T, C<->G transversions to A<->C, G<->T transversions
        const TypedDagNode< Simplex >*                      base_frequencies;
        
    };
    
}

#endif
