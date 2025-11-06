#include "ChronoToPhyloFunction.h"

#include <cstddef>

#include "RbBoolean.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ChronoToPhyloFunction::ChronoToPhyloFunction(const TypedDagNode<Tree> *c, const TypedDagNode<RbVector<double> > *r) : TypedFunction<Tree>( NULL ),
    chrono( c ),
    rates( r )
{
    // add the chrono parameter as a parent
    addParameter( chrono );
    addParameter( rates );
    
    update();
}


ChronoToPhyloFunction::~ChronoToPhyloFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ChronoToPhyloFunction* ChronoToPhyloFunction::clone( void ) const
{
    return new ChronoToPhyloFunction( *this );
}


// NB: Could be vastly improved with TreeListener events and touchSpecialization
void ChronoToPhyloFunction::update( void )
{
    delete value;
    
    // get a copy of the parent tree
    value = chrono->getValue().clone();
    
    // make the tree use branch lengths instead of ages
    value->getRoot().setUseAges(false, true);
    
    const RbVector<double> &this_rates = rates->getValue();
    
    size_t num_rates = this_rates.size();
    for (size_t i=0; i<num_rates; ++i)
    {
        double old_branch_length = chrono->getValue().getNode( i ).getBranchLength();
        double new_branch_length = old_branch_length * this_rates[i];
        value->getNode( i ).setBranchLength( new_branch_length );
    }
    
}


void ChronoToPhyloFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == chrono)
    {
        chrono = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if ( rates == oldP )
    {
        rates = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    
}

