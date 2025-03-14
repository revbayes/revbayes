//
//  FeatureInformedRateFunction.hpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/12/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#ifndef FeatureInformedRateFunction_hpp
#define FeatureInformedRateFunction_hpp


#include <vector>

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
    class DagNode;
    template <class valueType> class RbVector;
    template <class valueType> class TypedDagNode;
    
    class FeatureInformedRateFunction : public TypedFunction<RbVector<RbVector<double> > > {
        
    public:
        
        FeatureInformedRateFunction(
            const TypedDagNode< RbVector<RbVector<RbVector<long> > > >* cf,
            const TypedDagNode< RbVector<RbVector<RbVector<double> > > >* qf,
            const TypedDagNode< RbVector<double> >* cp,
            const TypedDagNode< RbVector<double> >* qp,
            const TypedDagNode< double >* nr
        );
        virtual                                                                         ~FeatureInformedRateFunction(void);
        
        // public member functions
        FeatureInformedRateFunction*                      clone(void) const;
        void                                              update(void);
        
    protected:
        
        void                                              swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        // members
        const TypedDagNode< RbVector<RbVector<RbVector<long> > > >*    categorical_features;
        const TypedDagNode< RbVector<RbVector<RbVector<double> > > >*  quantitative_features;
        const TypedDagNode< RbVector<double> >*  categorical_params;
        const TypedDagNode< RbVector<double> >*  quantitative_params;
        const TypedDagNode< double >* null_rate;
        
        size_t numCategoricalFeatures;
        size_t numQuantitativeFeatures;
        size_t numDim1;
        size_t numDim2;
        
    };
}


#endif /* FeatureInformedRateFunction_hpp */
