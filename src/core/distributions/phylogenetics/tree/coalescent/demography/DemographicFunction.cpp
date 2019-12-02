#include "DemographicFunction.h"
#include "DagNode.h"
#include "RbException.h"

using namespace RevBayesCore;

DemographicFunction::DemographicFunction(void) : Cloneable()
{
    
}


DemographicFunction::DemographicFunction(const DemographicFunction &f) :
    variables( f.variables )
{
    
    for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
    {
        const DagNode *the_node = *it;
        
        // add myself to the set of monitors
//        the_node->addMonitor( this );
        
        // tell the node that we have a reference to it (avoids deletion)
        the_node->incrementReferenceCount();
    }
    
}


DemographicFunction::~DemographicFunction(void)
{
    
    for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
    {
        const DagNode *the_node = *it;
        
        // remove myself to the set of monitors
//        the_node->removeMonitor( this );
        
        
        // tell the node that we have a reference to it (avoids deletion)
        if ( the_node->decrementReferenceCount() == 0 )
        {
            delete *it;
        }
    }
    
}


DemographicFunction& DemographicFunction::operator=(const DemographicFunction &f)
{
    
    if ( this != &f )
    {
        
        for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); ++it)
        {
            const DagNode *the_node = *it;
            
            // remove myself to the set of monitors
//            the_node->removeMonitor( this );
            
            // tell the node that we have a reference to it (avoids deletion)
            if ( the_node->decrementReferenceCount() == 0 )
            {
                delete *it;
            }
        }
        
        // set the nodes (we don't own them)
        variables = f.variables;
        
        for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); ++it)
        {
            
            const DagNode *the_node = *it;
            
            // add myself to the set of monitors
//            the_node->addMonitor( this );
            
            // tell the node that we have a reference to it (avoids deletion)
            the_node->incrementReferenceCount();
        }
    }
    
    return *this;
}


void DemographicFunction::addVariable(const DagNode *n)
{
    
    variables.push_back( n );
    
    // add myself to the set of monitors
//    n->addMonitor( this );
    
    // tell the node that we have a reference to it (avoids deletion)
    n->incrementReferenceCount();
    
    // now we need to sort again
//    sortNodesByName();
    
}


const std::vector<const DagNode *>& DemographicFunction::getDagNodes(void) const
{
    
    return variables;
}



//void DemographicFunction::removeVariable(DagNode *n)
//{
//
//    for (std::vector<DagNode*>::iterator it = variables.begin(); it != variables.end(); ++it)
//    {
//        if ( *it == n )
//        {
//            // move myself to the set of monitors
////            n->removeMonitor( this );
//
//            // tell the node that we have a reference to it (avoids deletion)
//            if ( n->decrementReferenceCount() == 0 )
//            {
//                delete *it;
//            }
//
//            variables.erase( it );
//
//            // we assume that nodes are unique
//            break;
//        }
//    }
//
//    // now we need to sort again
////    sortNodesByName();
//
//}


//void DemographicFunction::setDagNodes( const std::vector<DagNode *> &args)
//{
//
//    for (std::vector<DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
//    {
//        DagNode *the_node = *it;
//
//        // remove myself to the set of monitors
////        the_node->removeMonitor( this );
//
//        // tell the node that we have a reference to it (avoids deletion)
//        if ( the_node->decrementReferenceCount() == 0 )
//        {
//            delete *it;
//        }
//    }
//
//    // set the nodes (we don't own them)
//    variables = args;
//
//    for (std::vector<DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
//    {
//
//        DagNode *the_node = *it;
//
//        // add myself to the set of monitors
////        the_node->addMonitor( this );
//
//        // tell the node that we have a reference to it (avoids deletion)
//        the_node->incrementReferenceCount();
//    }
//
////    sortNodesByName();
//
//}


void DemographicFunction::swapNode(const DagNode *old_node, const DagNode *new_node)
{

    // error catching
    std::vector<const DagNode*>::iterator it = find(variables.begin(), variables.end(), old_node);
    
    if (it == variables.end())
    {
        throw RbException("Cannot replace DAG node with name\"" + old_node->getName() + "\" in this demographic function because the demographic function doesn't hold this DAG node.");
    }
    
    // remove myself from the old node and add myself to the new node
//    oldN->removeMonitor( this );
//    newN->addMonitor( this );
    
    // increment and decrement the reference counts
    new_node->incrementReferenceCount();
    if ( old_node->decrementReferenceCount() == 0 )
    {
        throw RbException("Memory leak in demographic function. Please report this bug to Sebastian.");
    }
    
    it = variables.insert( it, new_node );
    variables.erase( it + 1 );
    
    // now delegate to the derive classes
    swapNodeInternal(old_node, new_node);
}


std::ostream& operator<<(std::ostream& o, const DemographicFunction& x)
{
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const DemographicFunction& x)
{
    o << "DemographicFunction";
    
    return o;
}
