#include "ExtantTreeFunction.h"

#include <cassert>
#include <cstddef>
#include <vector>

#include "RbBoolean.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ExtantTreeFunction::ExtantTreeFunction(const TypedDagNode<Tree> *t) : TypedFunction<Tree>( new Tree() ),
    tau( t )
{
    // add the lambda parameter as a parent
    addParameter( tau );
    
    update();
}


ExtantTreeFunction::~ExtantTreeFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ExtantTreeFunction* ExtantTreeFunction::clone( void ) const
{
    return new ExtantTreeFunction( *this );
}


void ExtantTreeFunction::reInitialized( void )
{
    ; // *value = tau->getValue();
    
}

// NB: Could be vastly improved with TreeListener events and touchSpecialization
void ExtantTreeFunction::update( void )
{
    delete value;
    
    // get a copy of the parent tree
    value = tau->getValue().clone();
    
    // Extant the parent tree
    pruneFossilsRecursively( &value->getRoot() );
    
    value->setRoot( &value->getRoot(), true );
}

void ExtantTreeFunction::pruneFossilsRecursively( TopologyNode *node )
{
    

    // fossil node
    if ( node->isFossil() == true )
    {
        
        // get the parent and grand parent
        TopologyNode *parent = &node->getParent();
        TopologyNode *new_child = NULL;
        TopologyNode *new_parent = NULL;
        TopologyNode *recursive_node = NULL;

        // tip node
        if ( node->isTip() == true )
        {
            
            new_child = &parent->getChild(0);
            if ( new_child == node )
            {
                new_child = &parent->getChild(1);
            }

            // we must remove the child from its old parent
            parent->removeChild( new_child );
        }
        else
        {
            assert( node->getNumberOfChildren() == 1 );
            
            new_child = &node->getChild(0);
            recursive_node = new_child;

            // we must remove the child from its old parent
            node->removeChild( new_child );
        }
        
        // check if the parent is the root node
        if ( parent->isRoot() == true )
        {
            new_child->setParent( NULL );
            this->value->setRoot( new_child, true );
        }
        else
        {
            new_parent = &parent->getParent();
            new_parent->removeChild( parent );
            
            new_parent->addChild( new_child );
            new_child->setParent( new_parent );
        }
        

        if ( recursive_node != NULL )
        {
            pruneFossilsRecursively( recursive_node );
        }
        
    }
    else
    {
        
        std::vector<TopologyNode*> children = node->getChildren();
        for (size_t i = 0; i < children.size(); ++i)
        {
            pruneFossilsRecursively( children[i] );
        }
        
    }

    
}


void ExtantTreeFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}

