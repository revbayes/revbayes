#ifndef DispersalExtinctionRootStructureFunction_H
#define DispersalExtinctionRootStructureFunction_H

#include <cstddef>
#include <vector>
#include <map>

#include "Simplex.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class DispersalExtinctionRootStructureFunction : public TypedFunction<Simplex> {
        
    public:
        DispersalExtinctionRootStructureFunction(TypedDagNode<RbVector<double> >* rf, TypedDagNode<Simplex>* rs);
        virtual                                            ~DispersalExtinctionRootStructureFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        DispersalExtinctionRootStructureFunction*           clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(const DagNode* affecter);
        void                                                restore(const DagNode *restorer);
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(const DagNode *toucher );
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        void                                                makeBits(void);
        void                                                makeIdxByRangeSize(void);
        
        // members
        const TypedDagNode<RbVector<double> >*              root_frequencies;
        const TypedDagNode<Simplex>*                        rangeSize;
//        const TypedDagNode<std::int64_t>*                            maxRangeSize;
        
        std::vector<std::vector<unsigned> >                 bits;
        std::map<std::vector<unsigned>, unsigned>           inverseBits;
        std::vector<std::vector<unsigned> >                 idxByRangeSize;
        size_t                                              numCharacters;
        size_t                                              num_states;
    };
    
}


#endif /* defined(__revbayes_proj__DispersalExtinctionRootStructureFunction__) */
