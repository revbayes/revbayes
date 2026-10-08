#include "DeterministicNodeBase.h"

#include "DagNode.h"
#include "DynamicNodeBase.h"
#include "Function.h"
#include "RbException.h"

using namespace RevBayesCore;

/** Take ownership of a function and initialize its lazy-update state. */
DeterministicNodeBase::DeterministicNodeBase(Function *f) :
    function( f ),
    needs_update( true ),
    force_update( f->forceUpdates() )
{
}


/** Clone the owned function but leave graph attachment to the typed node. */
DeterministicNodeBase::DeterministicNodeBase(const DeterministicNodeBase &n) :
    function( n.function->clone() ),
    needs_update( true ),
    force_update( n.function->forceUpdates() )
{
}


/** Destroy the function owned by this implementation base. */
DeterministicNodeBase::~DeterministicNodeBase(void)
{
    delete function;
}


/** Replace the owned function and rebuild its incoming DAG edges for the existing node. */
void DeterministicNodeBase::assign(const DeterministicNodeBase &n, DagNode &owner)
{
    // Remove us as the child of the function parameters
    detachFromFunctionParameters( owner );

    // Delete the function
    delete function;

    // Recreate the function
    function = n.function->clone();

    // Get the parameters from the new function and add us as child of them in the DAG
    attachToFunctionParameters( owner );

    needs_update = true;
    force_update = function->forceUpdates();
}


/** Register the deterministic node as a child of every function parameter. */
void DeterministicNodeBase::attachToFunctionParameters(DagNode &owner) const
{
    // Get the parameters from the function and add us as a child of them in the DAG
    const std::vector<const DagNode*> &parents = function->getParameters();
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->addChild( &owner );

        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }
}


/** Remove the deterministic node's incoming edges and release its references to its parameters. */
void DeterministicNodeBase::detachFromFunctionParameters(DagNode &owner) const
{
    // Remove us as the child of the function parameters
    // Copy the vector because releasing a parameter may delete it and invalidate related storage.
    std::vector<const DagNode*> parents = function->getParameters();
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->removeChild( &owner );

        // Decrement the reference count and check whether we need to delete the DAG node
        if ( (*it)->decrementReferenceCount() == 0 )
        {
            delete *it;
        }
    }
}


/** A deterministic node passes affected-node discovery directly to its descendants. */
void DeterministicNodeBase::getAffected(DagNode &owner, RbOrderedSet<DagNode*> &affected) const
{
    owner.getAffectedNodes( affected );
}


/** Collect integrated stochastic ancestors by delegating through every current parent. */
void DeterministicNodeBase::getIntegratedParents(const DagNode &owner, RbOrderedSet<DagNode*> &integratedParents) const
{
    std::vector<const DagNode*> parents = owner.getParents();                                                                     //!< Get the set of parents (empty set here)

    // delegate up the DAG
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->getIntegratedParents( integratedParents );
    }
}


/**
 * Get the parents of this node. Simply ask the function to provide its parameters,
 * no need to keep parents here.
 */
std::vector<const DagNode*> DeterministicNodeBase::getParents(void) const
{
    return function->getParameters();
}


/** A deterministic node is constant exactly when every function parameter is constant. */
bool DeterministicNodeBase::isConstant(void) const
{
    // iterate over all parents and only if all parents are constant then this node is constant too
    const std::vector<const DagNode*> &parents = function->getParameters();
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        if ( (*it)->isConstant() == false )
        {
            return false;
        }
    }

    return true;
}


/**
 * Keep the current value of the node.
 * At this point, we just delegate to the children.
 */
void DeterministicNodeBase::keepMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *affecter)
{
    // delegate call to base class
    // this will unset the touched flag if it was set
    dynamicState.keepMe();

    // allow specialized recovery in functions
    function->keep( affecter );

    // delegate call
    owner.keepAffected();

    // clear the list of touched element indices
    owner.clearTouchedElementIndices();
}


/** Notify the function that its model has been reinitialized. */
void DeterministicNodeBase::reInitializeMe(void)
{
    function->reInitialized();
}


/** Restore the old value of the node and tell affected. */
void DeterministicNodeBase::restoreMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *restorer)
{
    // The value has been changed so we need to flag for recomputing the value. We need to do that even
    // if the touched flag is unset because a reset call from one parameter may have unset it before
    // another parameter was restored. We therefore update our value whenever a parameter is restored.
    needs_update = true;

    // we just mark ourselves as clean, albeit perhaps not being updated
    dynamicState.restoreMe();

    // call for potential specialized handling (e.g. internal flags)
    function->restore( restorer );

    // clear the list of touched element indices
    owner.clearTouchedElementIndices();

    // delegate call
    owner.restoreAffected();
}


/**
 * This function replaces the earlier swapParameter function. If we rely on the
 * internal RevBayesCore::Function to manage our parents, we simply need to ask
 * the function to swap its parameters, and then manage the connection of the
 * parents (parameters) to this node.
 */
void DeterministicNodeBase::swapParent(DagNode &owner, const DagNode *oldParent, const DagNode *newParent)
{
    // We are sure to get into trouble if either one of these is NULL
    if ( oldParent == NULL || newParent == NULL )
    {
        throw RbException( "Attempt to swap NULL function parameter of RevBayesCore::DeterministicNode" );
    }

    // This throws an error if the oldParent cannot be found
    function->swapParameter( oldParent, newParent );

    oldParent->removeChild( &owner );
    if ( oldParent->decrementReferenceCount() == 0 )
    {
        delete oldParent;
    }

    newParent->addChild( &owner );
    newParent->incrementReferenceCount();
    owner.touch();
}


/** Touch this node for recalculation. */
void DeterministicNodeBase::touchMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *toucher, bool touchAll)
{
    // store if the state of the variable was dirty (needed an update)
    bool needed_update = needs_update;
    bool was_touched = dynamicState.isTouched();

    // delegate call to base class; this will set the touched flag if it wasn't set already
    dynamicState.touchMe();

    // We need to touch the function always because of specialized touch functionality in some
    // functions, like vector functions. This is essential for lazy evaluation.
    function->touch( toucher );

    // mark for update
    needs_update = true;

    // The condition is historically always true; retaining it preserves the current propagation structure.
    if ( needed_update == false || was_touched == false || true )
    {
        // Dispatch the touch message to downstream nodes
        owner.touchAffected( touchAll );
    }
}
