//
//  DistanceRateModifierFunction.h
//  rb_mlandis
//
//  Created by Michael Landis on 5/17/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__DistanceRateModifierFunction__
#define __rb_mlandis__DistanceRateModifierFunction__

#include "TimeAtlas.h"
#include "CharacterHistoryRateModifier.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    class DistanceRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        DistanceRateModifierFunction(const TypedDagNode<double>* dp, const TimeAtlas* ta,  bool uadj=false, bool uav=false, bool udd=false);
        virtual ~DistanceRateModifierFunction(void);                                                                                                  //!< Virtual destructor
        
        // public member functions
        DistanceRateModifierFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                                    update(void);
        
    protected:
        void                                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                                             distancePower;
        
        
    };
}


#endif /* defined(__rb_mlandis__DistanceRateModifierFunction__) */
