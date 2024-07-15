/**
 * This file contains the declaration of the PoMo stationary frequencies funcion class.
 * This class is derived from the function class and is used to
 * compute the stationary frequencies of the PoMo model.
 */



#ifndef PoMoReversibleMutationRatesFunction_H
#define PoMoReversibleMutationRatesFunction_H

#include <vector>

#include "Simplex.h"
#include "TypedFunction.h"
#include "RateGenerator.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class PoMoReversibleMutationRatesFunction : public TypedFunction<  RbVector<double> > {
        
    public:
        PoMoReversibleMutationRatesFunction( long na, const TypedDagNode< Simplex > *bf, const TypedDagNode< RbVector<double> > *ex);
        

        virtual                                            ~PoMoReversibleMutationRatesFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoReversibleMutationRatesFunction*                clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        long                                                K;
        const TypedDagNode< Simplex > *                     pi;
        const TypedDagNode< RbVector<double> > *            rho;
        
    };
    
}

#endif
