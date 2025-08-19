#include "FeatureInformedRateFunction.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

#include <cmath>
#include <cstddef>

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


FeatureInformedRateFunction::FeatureInformedRateFunction(
    const TypedDagNode< RbVector<RbVector<RbVector<std::int64_t> > > >* cf,
    const TypedDagNode< RbVector<RbVector<RbVector<double> > > >* qf,
    const TypedDagNode< RbVector<double> >* cp,
    const TypedDagNode< RbVector<double> >* qp) :
        TypedFunction<RbVector<RbVector<double> > >( new RbVector<RbVector<double> >() ),
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
    
    numCategoricalFeatures = categorical_features->getValue().size();
    numQuantitativeFeatures = quantitative_features->getValue().size();
    numDim1 = categorical_features->getValue()[0].size();
    numDim2 = categorical_features->getValue()[0][0].size();
    
    // refresh values
//    this->setForceUpdates(true);
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
    
    // get all parent node values
    const RbVector<RbVector<RbVector<std::int64_t> > >& cf = categorical_features->getValue();
    const RbVector<RbVector<RbVector<double> > >& qf = quantitative_features->getValue();
    const RbVector<double>& cp = categorical_params->getValue();
    const RbVector<double>& qp = quantitative_params->getValue();

    // initialize new relative rates (=1)
    RbVector<RbVector<double> > rates(numDim1, RbVector<double>(numDim2, 1.0));
    
    // apply categorical scalers
    for (size_t i = 0; i < numCategoricalFeatures; i++) {
        double cp_i = std::exp(cp[i]);
        for (size_t j = 0; j < numDim1; j++) {
            for (size_t k = 0; k < numDim2; k++) {
                if (cf[i][j][k] != 0) {
                    rates[j][k] *= cp_i;
                }
            }
        }
    }
    // apply quantitative scalers
    for (size_t i = 0; i < numQuantitativeFeatures; i++) {
        for (size_t j = 0; j < numDim1; j++) {
            for (size_t k = 0; k < numDim2; k++) {
//                rates[j][k] *= std::pow( qf[i][j][k], qp[i] );
                if (std::isnan(qf[i][j][k]) == false) {
                    rates[j][k] *= std::exp( qf[i][j][k] * qp[i] );
                }
            }
        }
    }
    
    // set new values
    (*value) = rates;
}


void FeatureInformedRateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == categorical_features)
    {
        categorical_features = static_cast<const TypedDagNode< RbVector<RbVector<RbVector<std::int64_t> > > >* >( newP );
    }
    if (oldP == quantitative_features)
    {
        quantitative_features = static_cast<const TypedDagNode< RbVector<RbVector<RbVector<double> > > >* >( newP );
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
