#ifndef DispersalExtinctionRateStructureFunction_H
#define DispersalExtinctionRateStructureFunction_H

#include <cstddef>
#include <vector>
#include <map>

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class DispersalExtinctionRateStructureFunction : public TypedFunction<RbVector<double> > {
        
    public:
        DispersalExtinctionRateStructureFunction(const TypedDagNode<RbVector<RbVector<double> > >* dr, TypedDagNode<RbVector<double> >* er, TypedDagNode<std::int64_t>* rs);
        virtual                                            ~DispersalExtinctionRateStructureFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        DispersalExtinctionRateStructureFunction*           clone(void) const;                                                                  //!< Create an independent clone
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        void                                                makeBits(void);
        void                                                makeTransitions(void);
        
        // members
        const TypedDagNode<RbVector<RbVector<double> > >*   dispersalRates;
        const TypedDagNode<RbVector<double> >*              extinctionRates;
        const TypedDagNode<std::int64_t>*                            maxRangeSize;
        
        std::vector<std::vector<unsigned> >                 bits;
        std::map<std::vector<unsigned>, unsigned>           inverseBits;
        std::vector<std::vector<unsigned> >                 transitions;
        std::vector<std::vector<unsigned> >                 lossOrGain;
        std::vector<std::vector<std::vector<unsigned> > >   transitionAreas;
        size_t                                              numCharacters;
        size_t                                              num_states;
    };
    
}


#endif /* defined(__revbayes_proj__DispersalExtinctionRateStructureFunction__) */
