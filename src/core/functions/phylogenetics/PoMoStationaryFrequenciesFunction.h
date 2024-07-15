/**
 * This file contains the declaration of the PoMo stationary frequencies funcion class.
 * This class is derived from the function class and is used to
 * compute the stationary frequencies of the PoMo model.
 */



#ifndef PoMoStationaryFrequenciesFunction_H
#define PoMoStationaryFrequenciesFunction_H

#include <vector>

#include "Simplex.h"
#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class PoMoStationaryFrequenciesFunction : public TypedFunction< Simplex > {
        
    public:
        PoMoStationaryFrequenciesFunction( long na, long nv, const TypedDagNode< double > *ne, const TypedDagNode< Simplex > *bf, const TypedDagNode< RbVector<double> > *ex);
        

        virtual                                            ~PoMoStationaryFrequenciesFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoStationaryFrequenciesFunction*                  clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        long                                                K;
        long                                                V;
        const TypedDagNode< double > *                      N;
        const TypedDagNode< Simplex > *                     pi;
        const TypedDagNode< RbVector<double> > *            rho;
        
    };
    
}

#endif
