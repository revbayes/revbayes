//
//  FeatureInformedRateFunction.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/12/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "FeatureInformedRateFunction.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

#include <stddef.h>

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


//TypedFunction<MatrixReal>( new MatrixReal( mc + 1, (mc + 1) * (mc + 1), 0.0 ) ),
FeatureInformedRateFunction::FeatureInformedRateFunction(
    const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<long> > >* cf,
    const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* qf,
    const TypedDagNode< RevBayesCore::RbVector<double> >* cp,
    const TypedDagNode< RevBayesCore::RbVector<double> >* qp) :
        TypedFunction<RbVector<double> >( new RbVector<double>() ),
        categorical_features(cf),
        quantitative_features(qf),
        categorical_params(cp),
        quantitative_params(qp)

{
    // add parameters
    addParameter( categorical_features );
    addParameter( quantitative_features );
    addParameter( categorical_params );
    addParameter( quantitative_params );
    
    numStates = categorical_features->getValue()[0].size();
    
    // refresh values
    update();
}


FeatureInformedRateFunction::~FeatureInformedRateFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


FeatureInformedRateFunction* FeatureInformedRateFunction::clone( void ) const
{
    return new FeatureInformedRateFunction( *this );
}

void FeatureInformedRateFunction::update( void )
{
    RbVector<double> rates(numStates, 1.0);
    
    const RbVector<RbVector<long> >& cf = categorical_features->getValue();
    const RbVector<double>& cp = categorical_params->getValue();
    for (size_t i = 0; i < cf.size(); i++) {
        double cp_i = std::exp(cp[i]);
        for (size_t j = 0; j < cf[i].size(); j++) {
            if (cf[i][j] != 0) {
                rates[j] *= cp_i;
            }
        }
    }
    const RbVector<RbVector<double> >& qf = quantitative_features->getValue();
    const RbVector<double>& qp = quantitative_params->getValue();
    for (size_t i = 0; i < qf.size(); i++) {
        
        for (size_t j = 0; j < qf[i].size(); j++) {
            rates[j] *= std::pow( qf[i][j], qp[i] );
        }
    }
    
    (*value) = rates;
}


void FeatureInformedRateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == categorical_features)
    {
        categorical_features = static_cast<const TypedDagNode< RbVector<RbVector<long> > >* >( newP );
    }
    if (oldP == quantitative_features)
    {
        quantitative_features = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    if (oldP == categorical_params)
    {
        categorical_params = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == quantitative_params)
    {
        quantitative_params = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    
}
