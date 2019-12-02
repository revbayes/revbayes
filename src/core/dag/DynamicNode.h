#ifndef DynamicNode_H
#define DynamicNode_H

#include "TypedDagNode.h"
#include <set>

namespace RevBayesCore {
    
    template<class valueType>
    class DynamicNode : public TypedDagNode<valueType> {
        
    public:
        
        DynamicNode(const std::string &n);
        DynamicNode(const DynamicNode &n);
        virtual                                            ~DynamicNode(void);                                                              //!< Virtual destructor
        
        // pure virtual methods
        virtual DynamicNode<valueType>*                     clone(void) const = 0;
        
        // public methods
        virtual DagNode*                                    cloneDAG(DagNodeMap &nodesMap, std::map<std::string, const DagNode* > &names) const; //!< Clone the entire DAG which is connected to this node
        
        // this function provided for derived classes used in the language layer, which need to override it
        virtual const std::string&                          getRevTypeOfValue(void);                                                        //!< Get Rev language type of value
        
    protected:
        virtual void                                        keepMe(DagNode* affecter);                                                      //!< Keep value of this and affected nodes
        virtual void                                        restoreMe(DagNode *restorer);                                                   //!< Restore value of this node
        virtual void                                        touchMe(DagNode *toucher, bool touchAll);                                       //!< Tell affected nodes value is reset
        
        
        // members
        bool                                                touched;
        
    };
    
}

#include "RbException.h"
#include "RbOptions.h"

template<class valueType>
RevBayesCore::DynamicNode<valueType>::DynamicNode( const std::string &n ) : TypedDagNode<valueType>( n ),
    touched( true )
{
    // nothing to do here
}


template<class valueType>
RevBayesCore::DynamicNode<valueType>::DynamicNode( const DynamicNode<valueType> &n ) : TypedDagNode<valueType>( n ),
    touched( true )
{
    // nothing to do here
}


template<class valueType>
RevBayesCore::DynamicNode<valueType>::~DynamicNode( void )
{
    // We don't own the parents and hence we don't delete them. The owner of the graph needs to do that. Just tell them we are gone forever ... :(
}


/** Clone the entire graph: clone children, swap parents */
template<class valueType>
RevBayesCore::DagNode* RevBayesCore::DynamicNode<valueType>::cloneDAG( DagNodeMap& newNodes, std::map<std::string, const DagNode* > &names ) const
{
    
    // Return our clone if we have already been cloned
    if ( newNodes.find( this ) != newNodes.end() )
    {
        return ( newNodes[ this ] );
    }
    
    // just for self checking purposes we keep track of the names for the variables we already cloned
    if ( this->name != "" )
    {
        // check if we already added a variable with this name
        std::map<std::string, const DagNode* >::const_iterator n = names.find( this->name );
        if ( n == names.end() )
        {
            // no, we haven't cloned a variable with this name before
            names[ this->name ] = this;
        }
        else
        {
#ifdef DEBUG_SEBASTIAN
            const DagNode *orgCopy = n->second;
            std::cerr << "Ptr to org:\t" << orgCopy << "\t\t --- \t\t Ptr to desc:\t" << this << std::endl;
            std::cerr << "Cloning a DAG node with name '" << this->name << "' again, doh! Please tell this to Sebastian because it's most likely a bug." << std::endl;
#endif
        }
    }
    
    // Get a shallow copy
    DynamicNode* copy = clone();
    
    // Add this node and its copy to the map
    newNodes[ this ] = copy;
    
    // Parent management is delegated to derived classes, so get the parents through their getParents function
    std::vector<const DagNode*> my_parents = this->getParents();
    
    // We need to remove the copy as a child of our parents in order to stop recursive calls to
    // cloneDAG on our copy, its copy, etc, when we call cloneDAG on our parents
    for ( std::vector<const DagNode*>::const_iterator i = my_parents.begin(); i != my_parents.end(); ++i )
    {
        const DagNode *theParam = (*i);
        
        theParam->removeChild( copy );
        
        theParam->decrementReferenceCount();
    }
    
    // Now replace the parents of the copy (which are now the same as our parents) with the parent clones
    for ( std::vector<const DagNode*>::const_iterator i = my_parents.begin(); i != my_parents.end(); ++i )
    {
        // Get the i-th parent
        const DagNode *theParam = (*i);
        
        // Get its clone. If we already have cloned this parent (parameter), then we will get the previously created clone
        DagNode* theParamClone = theParam->cloneDAG( newNodes, names );
        
        // Add the copy back as a child of this parent so that the swapping works
        theParam->addChild( copy );
        theParam->incrementReferenceCount();
        
        // Swap the parent of the copy with its clone. This will remove the copy again as the child of our parent.
        copy->swapParent( theParam, theParamClone);
    }
    
    // Make sure the children clone themselves
    std::vector<DagNode*> children_to_clone = this->getChildren();
    for ( std::vector<DagNode* >::const_iterator i = children_to_clone.begin(); i != children_to_clone.end(); i++ )
    {
        DagNode *the_node = *i;
        std::string n = the_node->getName();
        the_node->cloneDAG( newNodes, names );
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
template<class valueType>
const std::string& RevBayesCore::DynamicNode<valueType>::getRevTypeOfValue(void)
{
    throw RbException( "Rev language type of dynamic DAG node value not known" );
}


/**
 * Keep the current value of the node.
 * At this point, we also need to make sure we update the stored ln probability.
 */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::keepMe( DagNode* /*affecter*/ )
{
    
    // unset the touched flag
    touched = false;
    
}


/** 
 * Restore the old value of the node and tell affected 
 */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::restoreMe(DagNode * /*restorer*/)
{

    // unset the touched flag
    touched = false;
    
}

/** 
 * Touch this node for recalculation 
 */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::touchMe( DagNode * /*toucher*/, bool /*touchAll*/ )
{

    // set the touched flag
    touched = true;
    
}



#endif

