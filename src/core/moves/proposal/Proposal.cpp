#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "DagNode.h"
#include "Move.h"
#include "Proposal.h"
#include "RbException.h"


using namespace RevBayesCore;



Proposal::Proposal(double p) :
    nodes(),
    move( nullptr ),
    targetAcceptanceRate(p)
{
    
}


Proposal::Proposal(const Proposal &p)  :
    nodes( p.nodes ),
    move( nullptr ),
    targetAcceptanceRate( p.targetAcceptanceRate )
{
    
    for (std::vector<DagNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        (*it)->incrementReferenceCount();
    }
    
}


Proposal::~Proposal( void )
{
    
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        DagNode *the_node = nodes[i];
        if ( the_node->decrementReferenceCount() == 0 )
        {
            delete the_node;
        }
        
    }
    
}



Proposal& Proposal::operator=(const Proposal &p)
{
    
    if ( this != &p )
    {
        
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            DagNode *the_node = nodes[i];
            if ( the_node->decrementReferenceCount() == 0 )
            {
                delete the_node;
            }
            
        }
        nodes.clear();
        
        nodes = p.nodes;
        move = p.move;
        
        for (std::vector<DagNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
        {
            (*it)->incrementReferenceCount();
        }
        
    }
    
    return *this;
}



/**
 * Add this node to our set of nodes.
 */
void Proposal::addNode( DagNode *n )
{
    bool exists = false;
    
    for (std::vector<DagNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        if ( *it == n )
        {
            exists = true;
            break;
        }
    }
    
    // only add the node if it doesn't exist already
    if ( n != nullptr && exists == false )
    {
        if (n->isClamped() && this->allowClamped() == false)
        {
            std::cout << "! In class defined as: " << typeid(*this).name() << std::endl;
            throw RbException("Cannot add the clamped node '" + n->getName() + "' to proposal" + 
            (this->getProposalName() == "" ? "." : (" " + this->getLongProposalName())));
        }
        nodes.push_back( n );
    
        // increment reference count
        n->incrementReferenceCount();
    }
    
    // delegate to the move
    if ( move != nullptr )
    {
        move->addNode( n );
    }
    
}


const Move* Proposal::getMove( void ) const
{
    
    return move;
}



/**
 * Get the vector of nodes on which this proposal is working on.
 *
 * \return  Const reference to a vector of nodes pointer on which the proposal operates.
 */
const std::vector<DagNode*>& Proposal::getNodes( void ) const
{
    
    return nodes;
}


std::vector<DagNode*> Proposal::identifyNodesToTouch(void)
{

    return nodes;
}



/**
 * Remove this node from our list of nodes.
 */
void Proposal::removeNode( RevBayesCore::DagNode *n )
{
    for (std::vector<DagNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        if ( *it == n )
        {
            nodes.erase( it );
            break;
        }
    }
    
    // delegate to the move
    if ( move != NULL )
    {
        move->removeNode( n );
    }

    if ( n->decrementReferenceCount() == 0 )
    {
        delete n;
    }

}


void Proposal::setMove(Move *m)
{
    
    move = m;
    
}

/**
 * Swap the old node with a new one.
 * This will be called for example when the entire model graph is cloned or
 * when we replace a variable with the same name (re-assignment).
 * Here we update our set and delegate to the derived class.
 */
void Proposal::swapNode(DagNode *oldP, DagNode *newP)
{
    bool exists = false;
    
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        
        if ( nodes[i] == oldP )
        {
            nodes[i] = newP;
            exists = true;
            break;
        }
    }

    if ( exists == true )
    {
        swapNodeInternal( oldP, newP );
        
        // increment and decrement the reference counts
        newP->incrementReferenceCount();
        if ( oldP->decrementReferenceCount() == 0 )
        {
            throw RbException("Memory leak in Function. Please report this bug to Sebastian.");
        }
    }
    else
    {
        throw RbException() << "Could not find the Proposal parameter to be swapped: " << oldP->getName();
    }
    
}

std::string Proposal::getLongProposalName() const
{
    std::vector<std::string> node_names;
    for(auto node: nodes)
	node_names.push_back(node->getName());

    return getProposalName() + "(" + StringUtilities::join(node_names,",") +")";
}
