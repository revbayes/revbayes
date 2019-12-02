//
//  SampledCladogenesisRootFrequenciesFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 8/11/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef SampledCladogenesisRootFrequenciesFunction_h
#define SampledCladogenesisRootFrequenciesFunction_h


#include "MatrixReal.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore {
    
class SampledCladogenesisRootFrequenciesFunction : public TypedFunction< Simplex > {

        
    public:
        SampledCladogenesisRootFrequenciesFunction(const TypedDagNode< RateGenerator > *rm, const TypedDagNode< MatrixReal > *cp, const TypedDagNode< Tree > *t, const TypedDagNode<double>* r);
        virtual                                             ~SampledCladogenesisRootFrequenciesFunction(void); //!< Virtual destructor
        
        // public member functions
        SampledCladogenesisRootFrequenciesFunction*                   clone(void) const; //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP); //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RateGenerator >*                rateMatrix;
        const TypedDagNode< MatrixReal >*                   cladogenesisProbabilities;
        const TypedDagNode< Tree >*                         tree;
        const TypedDagNode< double >*                       clock;
        
        size_t                                              numStates;
        double                                              initialLengthForStationarity;
        double                                              tunedLengthForStationarity;
        double                                              absorptionTolerance;
        double                                              stationarityDifferenceTolerance;
        double                                              stationarityPositiveTolerance;
        bool                                                isSampledSpeciationTree;
    };
    
}

#endif /* SampledCladogenesisRootFrequenciesFunction_h */
