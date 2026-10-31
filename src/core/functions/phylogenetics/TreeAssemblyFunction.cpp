#include "TreeAssemblyFunction.h"

#include <cstddef>
#include <vector>

#include "RbException.h"
#include "DagNode.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TreeChangeEventHandler.h"
#include "TreeChangeEventMessage.h"
#include <set>
#include "TypedDagNode.h"


/* NOTE: keep/restore/update functions are very fragile.
 *
 * There have been a lot of MCMC problems where old values are not properly restored
 * after a rejected move:
 * - See https://github.com/revbayes/revbayes/pull/549
 *
 * At the same time, we want to avoid marking branch lengths are being changed
 * when they have not been changed, because this causes too much recalculation
 * of transition matrices and partial likelihoods.
 *
 * Complications:
 * 1. we create a new copy of the tree, and do not modify the input tree.
 * 2. restore and update can be called multiple times after a single MCMC move.
 * 3. restore and update might clear the list of changed branches (in touchedIndices).
 * 4. the branch-length vector might not be up-to-date with TreeAssemblyFunction::update() or ::restore() is called.
 * 5. the topology might change before update is called.
 * 6. some topology sources (dnEmpiricalTree / mvEmpiricalTree) delete the Tree and allocate a
 *    new one, which fires no tree change event, so a touch without one means every branch moved.
 *    - See https://github.com/revbayes/revbayes/issues/547
 */

using namespace RevBayesCore;

TreeAssemblyFunction::TreeAssemblyFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *b) : TypedFunction<Tree>( NULL ),
    tau( t ),
    brlen( b )
{
    if ( tau->getValue().getNumberOfNodes() - 1 != brlen->getValue().size() )
    {
        throw RbException() << "Number of branches (" << (tau->getValue().getNumberOfNodes() - 1) << ") does not match the number of branch lengths (" << brlen->getValue().size() << ")";
    }

    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( brlen );

    value = tau->getValue().clone();

    listenTo( &tau->getValue() );
    brlenFlagDirty = true;
    update();
}


TreeAssemblyFunction::TreeAssemblyFunction(const TreeAssemblyFunction &f) : TypedFunction<Tree>( f ),
    tau( f.tau ),
    brlen( f.brlen )
{

    // the base class has created a new value instance
    // we need to delete it here to avoid memory leaks
    delete value;

    value = tau->getValue().clone();

    listenTo( &tau->getValue() );
    brlenFlagDirty = true;
    update();
}


TreeAssemblyFunction::~TreeAssemblyFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

    listenTo( NULL );
    // we own the tree, so we can let the base class destructor delete it
}


/** Subscribe to the input tree, which dnEmpiricalTree can replace under us. */
void TreeAssemblyFunction::listenTo( const Tree *t )
{
    if ( heard_from == t ) return;
    if ( heard_from != NULL ) const_cast<Tree*>( heard_from )->getTreeChangeEventHandler().removeListener( this );
    heard_from = t;
    if ( heard_from != NULL ) const_cast<Tree*>( heard_from )->getTreeChangeEventHandler().addListener( this );
    tau_touched = true;
}


void TreeAssemblyFunction::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    heard_event = true;

    // an age or a length of the input's own, which our copy does not take from it
    if ( m == TreeChangeEventMessage::BRANCH_LENGTH ) return;

    // passed on by update(), once our copy has the new shape
    shape_moved = true;
    event_nodes.push_back( n.getIndex() );
}



TreeAssemblyFunction* TreeAssemblyFunction::clone( void ) const
{
    return new TreeAssemblyFunction( *this );
}


void TreeAssemblyFunction::keep( const DagNode *affecter )
{
    //delegate to base class
    TypedFunction< Tree >::keep( affecter );

    // SH (20200221): This needs to stay, otherwise the MCMC does not update properly.
    touchedNodeIndices.clear();
    brlenFlagDirty = false;
    // shape_moved remembers what this drops, since update() may not have read it yet
    event_nodes.clear();

    // SH (20200221): We must not call update because it breaks the MCMC!
    // SH (20190913): There seems to be an issue if we use two replicates
    // So we need to make sure keep is only called after update!
//    update();

}


void TreeAssemblyFunction::reInitialized( void )
{

}


void TreeAssemblyFunction::restore( const DagNode *restorer )
{
    //delegate to base class
    TypedFunction< Tree >::restore( restorer );

    touchedNodeIndices.clear();
    brlenFlagDirty = false;

    // The flags are cleared first so this doesn't notify any listeners. They will have already
    // been notified by an upstream touch.
    update();
}


void TreeAssemblyFunction::touch(const DagNode *toucher)
{

    //delegate to base class
    TypedFunction< Tree >::touch( toucher );

    //reset flag
    brlenFlagDirty = true;

    if ( toucher == tau )
    {
        // drop stale per-branch indices so update() rewrites every branch length onto it
        touchedNodeIndices.clear();
        tau_touched = true;
    }

    if ( toucher == brlen )
    {
        const std::set<size_t> &touchedIndices = toucher->getTouchedElementIndices();
        touchedNodeIndices.insert(touchedIndices.begin(), touchedIndices.end());
    }

    if ( toucher != tau && touchedNodeIndices.size() == 0 && brlen != NULL )
    {
        for (size_t i = 0; i < brlen->getValue().size(); ++i)
        {
            touchedNodeIndices.insert(i);
        }
    }
}


/** Copy the input's topology onto our own tree, as a tree of branch lengths rather than one of ages. */
void TreeAssemblyFunction::copyTopology( void )
{
    // the assignment below sends a tree listener event
    std::set<TreeChangeEventListener*> listeners = value->getTreeChangeEventHandler().getListeners();
    for (TreeChangeEventListener *l : listeners) value->getTreeChangeEventHandler().removeListener( l );

    const Tree &in = tau->getValue();
    double root_length = in.getRoot().getBranchLength();

    // Tree's assignment keeps the node indices of the input, which the branch lengths are indexed by
    *value = in;
    // the root has no branch
    value->getRoot().setBranchLength( root_length, false );

    for (TreeChangeEventListener *l : listeners) value->getTreeChangeEventHandler().addListener( l );
}


void TreeAssemblyFunction::update( void )
{
    const std::vector<double> &v = brlen->getValue();

    // bind to the current topology
    listenTo( &tau->getValue() );

    // silence from a tree that was touched means it was replaced wholesale
    bool replaced = ( tau_touched == true && heard_event == false );
    if ( shape_moved == true || replaced == true )
    {
        copyTopology();
        if ( value->getNumberOfNodes() - 1 != v.size() )
        {
            throw RbException() << "Number of branches (" << (value->getNumberOfNodes() - 1) << ") does not match the number of branch lengths (" << v.size() << ")";
        }
    }

    // quietly: what a reader has to be told is the input's tree change events, passed on below
    for (size_t i = 0; i < v.size(); ++i)
    {
        value->getNode(i).setBranchLength( v[i], false );
    }

    // restore() clears this, and then there is nothing to notify a listener about
    if ( brlenFlagDirty == true )
    {
        if ( replaced == true )
        {
            for (size_t i = 0; i < value->getNumberOfNodes(); ++i)
            {
                value->getTreeChangeEventHandler().fire( value->getNode(i), TreeChangeEventMessage::TOPOLOGY );
            }
        }
        else
        {
            for (size_t j = 0; j < event_nodes.size(); ++j)
            {
                if ( event_nodes[j] < value->getNumberOfNodes() )
                {
                    value->getTreeChangeEventHandler().fire( value->getNode( event_nodes[j] ), TreeChangeEventMessage::TOPOLOGY );
                }
            }
        }

        // a branch length move does not touch the tree, so only its touched elements say
        for (std::set<size_t>::iterator it = touchedNodeIndices.begin(); it != touchedNodeIndices.end(); ++it)
        {
            if ( *it < v.size() )
            {
                value->getTreeChangeEventHandler().fire( value->getNode( *it ), TreeChangeEventMessage::BRANCH_LENGTH );
            }
        }
    }

    event_nodes.clear();
    shape_moved = false;
    heard_event = false;
    tau_touched = false;
    touchedNodeIndices.clear();
    brlenFlagDirty = false;
}



void TreeAssemblyFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        listenTo( &tau->getValue() );
        // our own copy keeps its identity, and update() copies the new topology onto it
        touchedNodeIndices.clear();
    }
    else if (oldP == brlen)
    {
        brlen = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}
