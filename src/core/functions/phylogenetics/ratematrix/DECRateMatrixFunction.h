#ifndef DECRateMatrixFunction_H
#define DECRateMatrixFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class DECRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        DECRateMatrixFunction(const TypedDagNode< RbVector<RbVector<double> > > *dr, const TypedDagNode< RbVector<RbVector<double> > > *er, const TypedDagNode< Simplex > *rs, bool cs=true, bool ex=false, bool uc=false, size_t mrs=0);
        virtual                                             ~DECRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        DECRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        size_t computeNumStates(size_t numAreas, size_t maxRangeSize, bool orderedStates);
        
        // members
        const TypedDagNode< RbVector<RbVector<double> > >*  dispersalRates;
        const TypedDagNode< RbVector<RbVector<double> > >*  extirpationRates;
        const TypedDagNode< Simplex >*                      rangeSize;
    };
    
}

#endif /* defined(__revbayes_proj__DECRateMatrixFunction__) */
