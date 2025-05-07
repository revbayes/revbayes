//
//  StateCountRateModifierFunction.hpp
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/16/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef StateCountRateModifierFunction_hpp
#define StateCountRateModifierFunction_hpp

#include <cstddef>

#include "CharacterHistoryRateModifier.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class StateCountRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        StateCountRateModifierFunction(const TypedDagNode<RbVector<double> >* sf, size_t nc);
        StateCountRateModifierFunction(const StateCountRateModifierFunction& m);
        virtual ~StateCountRateModifierFunction(void);                                                                                                  //!< Virtual destructor
        
        // public member functions
        StateCountRateModifierFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                                update(void);
        
    protected:
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<double> >*                              stateFactors;
        
    };
}

#endif /* StateCountRateModifierFunction_hpp */
