//
//  DispersalExtinctionRateStructureFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__DispersalExtinctionRateStructureFunction__
#define __revbayes_proj__DispersalExtinctionRateStructureFunction__

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class DispersalExtinctionRateStructureFunction : public TypedFunction<RbVector<double> > {
        
    public:
        DispersalExtinctionRateStructureFunction(const TypedDagNode<RbVector<RbVector<double> > >* dr, TypedDagNode<RbVector<double> >* er, TypedDagNode<long>* rs);
        virtual                                            ~DispersalExtinctionRateStructureFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        DispersalExtinctionRateStructureFunction*           clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(DagNode *toucher );
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        void                                                makeBits(void);
        void                                                makeTransitions(void);
        
        // members
        const TypedDagNode<RbVector<RbVector<double> > >*   dispersalRates;
        const TypedDagNode<RbVector<double> >*              extinctionRates;
        const TypedDagNode<long>*                            maxRangeSize;
        
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
