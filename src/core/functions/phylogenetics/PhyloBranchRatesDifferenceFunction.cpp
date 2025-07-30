#include "PhyloBranchRatesDifferenceFunction.h"

#include <cmath>

#include "DeterministicNode.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PhyloBranchRatesDifferenceFunction::PhyloBranchRatesDifferenceFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *d ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    tau( t ),
    branch_rates( d )
{
    // add the tau parameter as a parent
    addParameter( tau );
    addParameter( branch_rates );
       
    update();
}


PhyloBranchRatesDifferenceFunction::~PhyloBranchRatesDifferenceFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PhyloBranchRatesDifferenceFunction* PhyloBranchRatesDifferenceFunction::clone( void ) const
{
    return new PhyloBranchRatesDifferenceFunction( *this );
}


void PhyloBranchRatesDifferenceFunction::keep( const DagNode *affecter )
{
    //delegate to base class
    TypedFunction< RbVector<double> >::keep( affecter );
    
    // clear the indices, otherwise the MCMC does not update properly.
    touched_node_indices.clear();
    
}


void PhyloBranchRatesDifferenceFunction::reInitialized( void )
{
    
}


void PhyloBranchRatesDifferenceFunction::restore( const DagNode *restorer )
{
    //delegate to base class
    TypedFunction< RbVector<double> >::restore( restorer );
    
//    touched_node_indices.clear();
}


void PhyloBranchRatesDifferenceFunction::touch(const DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< RbVector<double> >::touch( toucher );
    
    //reset flag
//    touchedTopology = false;
    
    if ( toucher == branch_rates )
    {
        const std::set<size_t> &touched_indices = toucher->getTouchedElementIndices();
        touched_node_indices.insert(touched_indices.begin(), touched_indices.end());
    }
//    else if (toucher == tau)
//    {
//        touchedTopology = true;
//    }
}


void PhyloBranchRatesDifferenceFunction::update( void )
{
    
    RbVector<double> &v = *value;
    const Tree &tree = tau->getValue();
    
    size_t num_branches = 2*tree.getNumberOfTips() - 2;
    
    if ( v.size() != num_branches )
    {
        v.resize( num_branches );
    }
    
    const TopologyNode &root = tree.getRoot();
    recursiveUpdate(root);
    
}

void PhyloBranchRatesDifferenceFunction::recursiveUpdate( const TopologyNode& node )
{
    
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        const TopologyNode& parent = node.getParent();
        if ( parent.isRoot() )
        {
            (*this->value)[index] = 0.0;
        }
        else
        {
            size_t parent_index = parent.getIndex();
            double parent_value = branch_rates->getValue()[parent_index];
            double node_value = branch_rates->getValue()[index];
            (*this->value)[index] = node_value - parent_value;
        }
        
        // touch the index of the current
//        dag_node->addTouchedElementIndex(index);
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    
    for (size_t i = 0; i < num_children; ++i)
    {
        recursiveUpdate(node.getChild(i));
    }
    
}

void PhyloBranchRatesDifferenceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == branch_rates)
    {
        branch_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}



