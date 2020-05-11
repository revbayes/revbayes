//
//  HostSwitchRateModifierFunction.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/24/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef HostSwitchRateModifierFunction_h
#define HostSwitchRateModifierFunction_h


#include "CharacterHistoryRateModifier.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class HostSwitchRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        HostSwitchRateModifierFunction(const TypedDagNode<Tree>* t, const TypedDagNode<RbVector<double> >* s);
        HostSwitchRateModifierFunction(const HostSwitchRateModifierFunction& m);
        virtual ~HostSwitchRateModifierFunction(void);                                                                                                  //!< Virtual destructor
        
        // public member functions
        HostSwitchRateModifierFunction*                         clone(void) const;                                                          //!< Create an independent clone
        void                                                              update(void);
        
    protected:
        void                                                              swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                                tau;
        const TypedDagNode<RbVector<double> >*                            scale;
        
    };
}

#endif /* HostSwitchRateModifierFunction_h */
