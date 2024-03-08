#include "PhyloBranchRatesFunction.h"

#include <cmath>

#include "DeterministicNode.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PhyloBranchRatesFunction::PhyloBranchRatesFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *d, bool l ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    tau( t ),
    node_state( d ),
    as_log( l )
{
    // add the tau parameter as a parent
    addParameter( tau );
    addParameter( node_state );
       
    update();
}


PhyloBranchRatesFunction::~PhyloBranchRatesFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PhyloBranchRatesFunction* PhyloBranchRatesFunction::clone( void ) const
{
    return new PhyloBranchRatesFunction( *this );
}


void PhyloBranchRatesFunction::keep( const DagNode *affecter )
{
    //delegate to base class
    TypedFunction< RbVector<double> >::keep( affecter );
    
    // clear the indices, otherwise the MCMC does not update properly.
    touched_node_indices.clear();
    
}


void PhyloBranchRatesFunction::reInitialized( void )
{
    
}


void PhyloBranchRatesFunction::restore( const DagNode *restorer )
{
    //delegate to base class
    TypedFunction< RbVector<double> >::restore( restorer );
    
//    touched_node_indices.clear();
}


void PhyloBranchRatesFunction::touch(const DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< RbVector<double> >::touch( toucher );
    
    //reset flag
//    touchedTopology = false;
    
    if ( toucher == node_state )
    {
        const std::set<size_t> &touched_indices = toucher->getTouchedElementIndices();
        touched_node_indices.insert(touched_indices.begin(), touched_indices.end());
    }
//    else if (toucher == tau)
//    {
//        touchedTopology = true;
//    }
}


void PhyloBranchRatesFunction::update( void )
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

void PhyloBranchRatesFunction::recursiveUpdate( const TopologyNode& node )
{
    
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        size_t parent_index = node.getParent().getIndex();
        double parent_value = node_state->getValue()[parent_index];
        double node_value = node_state->getValue()[index];
        if ( as_log == true )
        {
            (*this->value)[index] = exp( (parent_value + node_value) / 2.0 );
        }
        else
        {
            (*this->value)[index] = (parent_value + node_value) / 2.0;
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

void PhyloBranchRatesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == node_state)
    {
        node_state = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}



