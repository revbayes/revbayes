//
//  RangeEvolutionRateModifierFunction.hpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef RangeEvolutionRateModifierFunction_hpp
#define RangeEvolutionRateModifierFunction_hpp

#include "CharacterHistoryRateModifier.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    class RangeEvolutionRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        RangeEvolutionRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<RbVector<RbVector<double> > >* c, size_t nc);
        RangeEvolutionRateModifierFunction(const RangeEvolutionRateModifierFunction& m);
        virtual ~RangeEvolutionRateModifierFunction(void);                                                                                                  //!< Virtual destructor
        
        // public member functions
        RangeEvolutionRateModifierFunction*                               clone(void) const;                                                          //!< Create an independent clone
        void                                                              update(void);
        
    protected:
        void                                                              swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                                       gainFactor;
        const TypedDagNode<double>*                                       lossFactor;
        const TypedDagNode<RbVector<RbVector<double> > >*                 context_matrix;
        
    };
}

#endif /* RangeEvolutionRateModifierFunction_hpp */
