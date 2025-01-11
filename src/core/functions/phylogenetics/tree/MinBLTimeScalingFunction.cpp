#include "MinBLTimeScalingFunction.h"

#include <vector>

#include "RbVector.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

MinBLTimeScalingFunction::MinBLTimeScalingFunction(const TypedDagNode< Tree >* tr, const TypedDagNode< RbVector<Taxon> >* tx, const TypedDagNode<double>* mbl) : TypedFunction< Tree >( new Tree() ),
    treeToTimeScale( tr ),
    taxonVector( tx ),
    minimumBranchLength( mbl )

{
    addParameter( treeToTimeScale );
    addParameter( taxonVector );
    addParameter( minimumBranchLength );
    
    update();
}


MinBLTimeScalingFunction::~MinBLTimeScalingFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


MinBLTimeScalingFunction* MinBLTimeScalingFunction::clone( void ) const
{
    
    return new MinBLTimeScalingFunction( *this );
}


void MinBLTimeScalingFunction::update( void )
{
    delete value;
    
    // get a copy of the parent tree
    value = treeToTimeScale->getValue().clone();
    
    // apply the function
    TreeUtilities::minBLTimeScaling( *value, taxonVector->getValue(), minimumBranchLength->getValue() );
}


void MinBLTimeScalingFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == treeToTimeScale)
    {
        treeToTimeScale = static_cast<const TypedDagNode< Tree >* >( newP );
    }
    else if (oldP == taxonVector)
    {
        taxonVector = static_cast<const TypedDagNode< RbVector<Taxon> >* >( newP );
    }
    else if (oldP == minimumBranchLength)
    {
        minimumBranchLength = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
