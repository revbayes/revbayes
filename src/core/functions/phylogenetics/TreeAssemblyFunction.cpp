#include "TreeAssemblyFunction.h"

#include <cstddef>
#include <vector>

#include "RbConstants.h"
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
 * 3. a branch length vector computed as a whole names no changed elements.
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

    // a kept edit is never undone, so no restore can use this record
    rebuilt_valid = false;
    before_move_valid = false;

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

    brlenFlagDirty = false;

    // The flags are cleared first so this doesn't notify any listeners. They will have already
    // been notified by an upstream touch.
    restoring = true;
    update();
    restoring = false;

    // The length vector may not be restored yet. Give our copy the lengths the readers went back to, or a
    // later move that proposes the same lengths again would find nothing changed.
    if ( before_move_valid == true )
    {
        for (size_t i = 0; i < before_move.size() && i < value->getNumberOfNodes(); ++i) value->getNode(i).setBranchLength( before_move[i], false );
        before_move_valid = false;
    }
}


void TreeAssemblyFunction::touch(const DagNode *toucher)
{

    //delegate to base class
    TypedFunction< Tree >::touch( toucher );

    // the lengths readers hold before this move, which a rejection returns them to
    if ( before_move_valid == false )
    {
        before_move.resize( value->getNumberOfNodes() );
        for (size_t i = 0; i < before_move.size(); ++i) before_move[i] = value->getNode(i).getBranchLength();
        before_move_valid = true;
    }

    //reset flag
    brlenFlagDirty = true;

    if ( toucher == tau )
    {
        // An assembly nothing reads never runs update(), so it can go many moves without
        // rewiring, and the record would describe an older shape than the one being undone.
        rebuilt_valid = false;
        tau_touched = true;
    }
}


/**
 * Move our own nodes into the input's shape instead of letting Tree's assignment rebuild them.
 * Index i is the same branch in both trees, so they can be diffed directly.
 * False where the trees do not correspond node for node, so nothing is mutated before everything
 * this touches has been validated and the caller can still copy.
 */
bool TreeAssemblyFunction::syncShape( void )
{
    const Tree &in = tau->getValue();
    const std::vector<TopologyNode*> &theirs = in.getNodes();
    const std::vector<TopologyNode*> &ours = value->getNodes();
    size_t n = ours.size();
    if ( n != theirs.size() or n == 0 ) return false;

    // After a rejected move our copy differs from the restored input only where we rewired it.
    bool full = ( restoring == false or rebuilt_valid == false );
    size_t count = ( full == true ? n : rebuilt.size() );

    // Only mvIndependentTopology moves the root to another node. Bail and let the caller copy.
    size_t root_index = in.getRoot().getIndex();
    if ( root_index >= n or value->getRoot().getIndex() != root_index ) return false;
    if ( ours[root_index]->getIndex() != root_index ) return false;

    for (size_t k = 0; k < count; ++k)
    {
        size_t i = ( full == true ? k : rebuilt[k] );
        if ( i >= n or ours[i]->getIndex() != i or theirs[i]->getIndex() != i ) return false;
        const std::vector<TopologyNode*> &tc = theirs[i]->getChildren();
        for (size_t c = 0; c < tc.size(); ++c)
        {
            size_t ci = tc[c]->getIndex();
            if ( ci >= n or ours[ci]->getIndex() != ci ) return false;
        }
    }

    std::vector<size_t> next;
    for (size_t k = 0; k < count; ++k)
    {
        size_t i = ( full == true ? k : rebuilt[k] );
        const std::vector<TopologyNode*> &tc = theirs[i]->getChildren();
        const std::vector<TopologyNode*> &oc = ours[i]->getChildren();
        bool same = ( tc.size() == oc.size() );
        for (size_t c = 0; same == true and c < tc.size(); ++c)
        {
            same = ( tc[c]->getIndex() == oc[c]->getIndex() );
        }
        if ( same == true ) continue;

        // quietly: update() passes on the input's own events once the shape is in place
        while ( ours[i]->getNumberOfChildren() > 0 ) ours[i]->removeChild( &ours[i]->getChild(0) );
        for (size_t c = 0; c < tc.size(); ++c)
        {
            TopologyNode *child = ours[ tc[c]->getIndex() ];
            ours[i]->addChild( child );
            child->setParent( ours[i], false );
        }
        next.push_back(i);
    }

    if ( full == true )
    {
        rebuilt = next;
        rebuilt_valid = true;
    }

    return true;
}


/** Copy the input's topology to our tree. Rewires in place where it can, and copies otherwise. */
void TreeAssemblyFunction::copyTopology( bool replaced )
{
    // the assignment below sends a tree listener event
    std::set<TreeChangeEventListener*> listeners = value->getTreeChangeEventHandler().getListeners();
    for (TreeChangeEventListener *l : listeners) value->getTreeChangeEventHandler().removeListener( l );

    const Tree &in = tau->getValue();
    double root_length = in.getRoot().getBranchLength();

    // A replaced tree can differ arbitrarily, so there is nothing for syncShape to diff against.
    if ( replaced == true or syncShape() == false )
    {
        rebuilt_valid = false;
        *value = in;
    }

    // Our branch lengths come from the vector, so the input's ages do not describe our branches.
    // Dropping them keeps the Newick writer from recomputing every length as parent age minus age.
    value->getRoot().setUseAges( false, true );

    // the root has no branch
    value->getRoot().setBranchLength( root_length, false );

    for (TreeChangeEventListener *l : listeners) value->getTreeChangeEventHandler().addListener( l );
}


void TreeAssemblyFunction::update( void )
{
    const std::vector<double> &v = brlen->getValue();

    // bind to the current topology
    listenTo( &tau->getValue() );

    // the lengths a reader last saw, before a copy below can replace our nodes
    std::vector<double> seen( v.size(), RbConstants::Double::nan );
    for (size_t i = 0; i < v.size() && i < value->getNumberOfNodes(); ++i) seen[i] = value->getNode(i).getBranchLength();

    // silence from a tree that was touched means it was replaced wholesale
    bool replaced = ( tau_touched == true && heard_event == false );
    if ( shape_moved == true || replaced == true )
    {
        copyTopology( replaced );
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

        // only a length that changed, whichever parent's touch announced it
        for (size_t i = 0; i < v.size(); ++i)
        {
            if ( v[i] != seen[i] ) value->getTreeChangeEventHandler().fire( value->getNode(i), TreeChangeEventMessage::BRANCH_LENGTH );
        }
    }

    event_nodes.clear();
    shape_moved = false;
    heard_event = false;
    tau_touched = false;
    brlenFlagDirty = false;
}



void TreeAssemblyFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        listenTo( &tau->getValue() );
        // our own copy keeps its identity, and update() copies the new topology onto it
    }
    else if (oldP == brlen)
    {
        brlen = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}
