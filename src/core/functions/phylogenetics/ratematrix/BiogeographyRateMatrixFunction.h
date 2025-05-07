#ifndef BiogeographyRateMatrixFunction_H
#define BiogeographyRateMatrixFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class BiogeographyRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        BiogeographyRateMatrixFunction(const TypedDagNode< RbVector<RbVector<double> > > *dr, const TypedDagNode<RbVector<double> > *er, size_t mrs=0);
        virtual                                             ~BiogeographyRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        BiogeographyRateMatrixFunction*                     clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        size_t computeNumStates(size_t numAreas, size_t maxRangeSize, bool orderedStates);
        
        // members
        const TypedDagNode< RbVector<RbVector<double> > >*  dispersalRates;
        const TypedDagNode< RbVector<double> >*  extirpationRates;
        
    };
    
}

#endif /* defined(__revbayes_proj__BiogeographyRateMatrixFunction__) */
