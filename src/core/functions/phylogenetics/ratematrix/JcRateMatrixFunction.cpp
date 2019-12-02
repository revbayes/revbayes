#include "JcRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

JcRateMatrixFunction::JcRateMatrixFunction(size_t ns) : TypedFunction<RateGenerator>( new RateMatrix_JC(ns) ) {
    
    update();
}


JcRateMatrixFunction::~JcRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



JcRateMatrixFunction* JcRateMatrixFunction::clone( void ) const {
    return new JcRateMatrixFunction( *this );
}


void JcRateMatrixFunction::update( void ) {
    // nothing to do here
}



void JcRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    // nothing to do
}


