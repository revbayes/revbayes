#include "FeatureInformedRateFunction.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

#include <cmath>
#include <climits>
#include <cstddef>

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


FeatureInformedRateFunction::FeatureInformedRateFunction(
    const TypedDagNode< RbVector<RbVector<RbVector<long> > > >* cf,
    const TypedDagNode< RbVector<RbVector<RbVector<double> > > >* qf,
    const TypedDagNode< RbVector<double> >* cp,
    const TypedDagNode< RbVector<double> >* qp,
    const TypedDagNode< double >* nr) :
        TypedFunction<RbVector<RbVector<double> > >( new RbVector<RbVector<double> >() ),
        categorical_features(cf),
        quantitative_features(qf),
        categorical_params(cp),
        quantitative_params(qp),
        null_rate(nr)

{
    // add parameters
    addParameter( categorical_features );
    addParameter( quantitative_features );
    addParameter( categorical_params );
    addParameter( quantitative_params );
    addParameter( null_rate );
    
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
    const RbVector<RbVector<RbVector<long> > >& cf = categorical_features->getValue();
    const RbVector<RbVector<RbVector<double> > >& qf = quantitative_features->getValue();
    const RbVector<double>& cp = categorical_params->getValue();
    const RbVector<double>& qp = quantitative_params->getValue();

    // initialize new relative rates (=1)
    RbVector<RbVector<double> > rates(numDim1, RbVector<double>(numDim2, 1.0));
    
    // First, compute product of relative rate effects while
    // ignoring elements encoded for null rate (e.g. nan)
    
    // using prod(x_i) is probably numerically unstable
    // convert to exp(log(sum(x_i)))
    
    // apply categorical scalers
    for (size_t i = 0; i < numCategoricalFeatures; i++) {
        double cp_i = std::exp(cp[i]);
        for (size_t j = 0; j < numDim1; j++) {
            for (size_t k = 0; k < numDim2; k++) {
//                if (std::isnan(cf[i][j][k]) == false) {
                if (cf[i][j][k] == 1) {
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
    
    // Next, revisit features and assign null rate to any
    // element encoded with null feature value (e.g. nan)
    
    // apply categorical scalers
    for (size_t i = 0; i < numCategoricalFeatures; i++) {
        double cp_i = std::exp(cp[i]);
        for (size_t j = 0; j < numDim1; j++) {
            for (size_t k = 0; k < numDim2; k++) {
//                if (std::isnan(cf[i][j][k])) {
                if (cf[i][j][k] < 0) {
                    rates[j][k] = null_rate->getValue();
                }
            }
        }
    }
    // apply quantitative scalers
    for (size_t i = 0; i < numQuantitativeFeatures; i++) {
        for (size_t j = 0; j < numDim1; j++) {
            for (size_t k = 0; k < numDim2; k++) {
                if (std::isnan(qf[i][j][k])) {
                    rates[j][k] = null_rate->getValue();
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
        categorical_features = static_cast<const TypedDagNode< RbVector<RbVector<RbVector<long> > > >* >( newP );
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
