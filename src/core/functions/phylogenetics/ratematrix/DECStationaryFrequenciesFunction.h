//
//  DECStationaryFrequenciesFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/4/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef DECStationaryFrequenciesFunction_h
#define DECStationaryFrequenciesFunction_h

#include "RateMatrix_DECRateMatrix.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class DECStationaryFrequenciesFunction : public TypedFunction<RbVector<double> > {
        
    public:
        DECStationaryFrequenciesFunction(const TypedDagNode<RateGenerator>* q);
        
        virtual                                             ~DECStationaryFrequenciesFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        DECStationaryFrequenciesFunction*                   clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RateGenerator>*  rateMatrix;
    };
    
}


#endif /* DECStationaryFrequenciesFunction_h */
