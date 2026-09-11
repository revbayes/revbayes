#include "DynamicNodeBase.h"

#include "DagNode.h"
#include "DagNodeMap.h"
#include "RbException.h"

#include <iostream>

using namespace RevBayesCore;

/** Initialize a new dynamic node in the touched state. */
DynamicNodeBase::DynamicNodeBase(void) :
    touched( true )
{
}


/** Copies intentionally begin touched, matching the historical dynamic-node copy behavior. */
DynamicNodeBase::DynamicNodeBase(const DynamicNodeBase & /*n*/) :
    touched( true )
{
}


/** Clone the entire graph: clone children, swap parents. */
DagNode* DynamicNodeBase::cloneDAG(const DagNode &owner, DagNodeMap &new_nodes, std::map<std::string, const DagNode*> &names) const
{
    // Return our clone if we have already been cloned
    if ( new_nodes.find( &owner ) != new_nodes.end() )
    {
        return new_nodes[ &owner ];
    }

    // just for self checking purposes we keep track of the names for the variables we already cloned
    if ( owner.getName() != "" )
    {
        // check if we already added a variable with this name
        std::map<std::string, const DagNode* >::const_iterator n = names.find( owner.getName() );
        if ( n == names.end() )
        {
            // no, we haven't cloned a variable with this name before
            names[ owner.getName() ] = &owner;
        }
        else
        {
#ifdef DEBUG_SEBASTIAN
            const DagNode *orgCopy = n->second;
            std::cerr << "Ptr to org:\t" << orgCopy << "\t\t --- \t\t Ptr to desc:\t" << &owner << std::endl;
            std::cerr << "Cloning a DAG node with name '" << owner.getName() << "' again, doh! Please tell this to Sebastian because it's most likely a bug." << std::endl;
#endif
        }
    }

    // Get a shallow copy
    DagNode *copy = owner.clone();

    // Add this node and its copy to the map
    new_nodes[ &owner ] = copy;

    // Parent management is delegated to derived classes, so get the parents through their getParents function
    std::vector<const DagNode*> my_parents = owner.getParents();

    // We need to remove the copy as a child of our parents in order to stop recursive calls to
    // cloneDAG on our copy, its copy, etc, when we call cloneDAG on our parents
    for ( std::vector<const DagNode*>::const_iterator i = my_parents.begin(); i != my_parents.end(); ++i )
    {
        const DagNode *the_param = *i;
        the_param->removeChild( copy );
        the_param->decrementReferenceCount();
    }

    // Now replace the parents of the copy (which are now the same as our parents) with the parent clones
    for ( std::vector<const DagNode*>::const_iterator i = my_parents.begin(); i != my_parents.end(); ++i )
    {
        const DagNode *the_param = *i;

        // Get its clone. If we already have cloned this parent (parameter), then we will get the previously created clone
        DagNode *the_param_clone = the_param->cloneDAG( new_nodes, names );

        // Add the copy back as a child of this parent so that the swapping works
        the_param->addChild( copy );
        the_param->incrementReferenceCount();

        // Swap the parent of the copy with its clone. This will remove the copy again as the child of our parent.
        copy->swapParent( the_param, the_param_clone );
    }

    // Make sure the children clone themselves
    std::vector<DagNode*> children_to_clone = owner.getChildren();
    for ( std::vector<DagNode*>::const_iterator i = children_to_clone.begin(); i != children_to_clone.end(); ++i )
    {
        (*i)->cloneDAG( new_nodes, names );
    }

    return copy;
}


/**
 * This function returns the Rev language type of the value. When used in the Rev
 * language layer, a DAG node must know the Rev language type of its value, or
 * construction of dynamic variables will not be safe. Here we just throw an
 * error, as a core DAG node need not know and should not know the language type
 * of its value.
 */
const std::string& DynamicNodeBase::getRevTypeOfValue(void) const
{
    throw RbException( "Rev language type of dynamic DAG node value not known" );
}


/** Report whether the dynamic node currently requires processing. */
bool DynamicNodeBase::isTouched(void) const
{
    return touched;
}


/**
 * Keep the current value of the node.
 * At this point, we also need to make sure we update the stored ln probability.
 */
void DynamicNodeBase::keepMe(void)
{
    // unset the touched flag
    touched = false;
}


/** Restore the old value of the node and tell affected. */
void DynamicNodeBase::restoreMe(void)
{
    // unset the touched flag
    touched = false;
}


/** Touch this node for recalculation. */
void DynamicNodeBase::touchMe(void)
{
    // set the touched flag
    touched = true;
}
