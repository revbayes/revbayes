#include "DemographicFunction.h"

#include <algorithm>
#include <string>

#include "DagNode.h"
#include "RbException.h"

using namespace RevBayesCore;

DemographicFunction::DemographicFunction(void) : Cloneable()
{
    
}

/**
 * @param[in]   f    The demographic function to copy.
 */
DemographicFunction::DemographicFunction(const DemographicFunction &f) :
    variables( f.variables )
{
    
    for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
    {
        const DagNode *the_node = *it;
        
        // tell the node that we have a reference to it (avoids deletion)
        the_node->incrementReferenceCount();
    }
    
}


DemographicFunction::~DemographicFunction(void)
{
    
    for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); it++)
    {
        const DagNode *the_node = *it; 
        
        // delete the node if there are no references to it
        if ( the_node->decrementReferenceCount() == 0 )
        {
            delete *it;
        }
    }
    
}

/**
 * @param[in]   f    The demographic function to copy.
 */
DemographicFunction& DemographicFunction::operator=(const DemographicFunction &f)
{
    
    if ( this != &f )
    {
        
        for (std::vector<const DagNode*>::iterator it = variables.begin(); it != variables.end(); ++it)
        {
            const DagNode *the_node = *it;
            
            // delete the node if there are no references to it
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
            
            // tell the node that we have a reference to it (avoids deletion)
            the_node->incrementReferenceCount();
        }
    }
    
    return *this;
}

/**
 * @param[in]   n    Pointer to the DAG node to be added to the nodes vector.
 */
void DemographicFunction::addVariable(const DagNode *n)
{
    
    variables.push_back( n );
      
    // tell the node that we have a reference to it (avoids deletion)
    n->incrementReferenceCount();   
}


/**
 * @return   The nodes vector containing the variables.
 */
const std::vector<const DagNode *>& DemographicFunction::getDagNodes(void) const
{
    
    return variables;
}


/**
 * @param[in]   old_node    Pointer to the DAG node to be replaced
 * @param[in]   new_node    Pointer to the DAG node replacing the other
 *
 */
void DemographicFunction::swapNode(const DagNode *old_node, const DagNode *new_node)
{

    // error catching
    std::vector<const DagNode*>::iterator it = find(variables.begin(), variables.end(), old_node);
    
    if (it == variables.end())
    {
        throw RbException() << "Cannot replace DAG node with name\"" <<  old_node->getName() << "\" in this demographic function because the demographic function doesn't hold this DAG node.";
    }
        
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
