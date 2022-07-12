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
    
    class FeatureInformedRateFunction : public TypedFunction<RbVector<double> > {
        
    public:
        
        FeatureInformedRateFunction(
            const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<long> > >* cf,
            const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* qf,
            const TypedDagNode< RevBayesCore::RbVector<double> >* cp,
            const TypedDagNode< RevBayesCore::RbVector<double> >* qp
        );
        virtual                                                                         ~FeatureInformedRateFunction(void);
        
        // public member functions
        FeatureInformedRateFunction*                      clone(void) const;
        void                                              update(void);
        
    protected:
        
        void                                              swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        // members
        const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<long> > >*    categorical_features;
        const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >*  quantitative_features;
        const TypedDagNode< RevBayesCore::RbVector<double> >*  categorical_params;
        const TypedDagNode< RevBayesCore::RbVector<double> >*  quantitative_params;
        size_t numStates;
    };
    
}


#endif /* FeatureInformedRateFunction_hpp */
