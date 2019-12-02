#include "InfiniteSitesRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

InfiniteSitesRateMatrixFunction::InfiniteSitesRateMatrixFunction(size_t ns) : TypedFunction<RateGenerator>( new RateMatrix_InfiniteSites(ns) )
{
    
    update();
}


InfiniteSitesRateMatrixFunction::~InfiniteSitesRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



InfiniteSitesRateMatrixFunction* InfiniteSitesRateMatrixFunction::clone( void ) const
{
    return new InfiniteSitesRateMatrixFunction( *this );
}


void InfiniteSitesRateMatrixFunction::update( void )
{
    // nothing to do here
}



void InfiniteSitesRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    // nothing to do
}


