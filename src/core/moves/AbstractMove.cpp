#include "AbstractMove.h"

#include <iterator>
#include <cstdint>
#include <string>

#include "DagNode.h"
#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;


/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    weight   The weight how often the proposal will be used (per iteration).
 * \param[in]    tuning   If auto tuning should be used.
 */
AbstractMove::AbstractMove( double weight, bool tuning ) :
    nodes(  ),
    affected_nodes(  ),
    weight( weight ),
    auto_tuning( tuning ),
    num_tried_current_period( 0 ),
    num_tried_total( 0 )
{
    
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    nodes    The weight how often the proposal will be used (per iteration).
 * \param[in]    weight   The weight how often the proposal will be used (per iteration).
 * \param[in]    tuning   If auto tuning should be used.
 */
AbstractMove::AbstractMove( const std::vector<DagNode*> &nodes, double weight, bool tuning ) :
    nodes( nodes ),
    affected_nodes( ),
    weight( weight ),
    auto_tuning( tuning ),
    num_tried_current_period( 0 ),
    num_tried_total( 0 )
{
    
    for (std::vector<DagNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        // get the pointer to the current node
        DagNode* the_node = *it;
        
        // add myself to the set of moves
        the_node->addMove( this );
        
        // increase the DAG node reference count because we also have a pointer to it
        the_node->incrementReferenceCount();
        
    }
    
    
    // remove all "core" nodes from affectedNodes so their probabilities are not double-counted
    for (size_t i = 0; i < affected_nodes.size(); ++i)
    {
        RbOrderedSet<DagNode*>::iterator it = affected_nodes.begin();
        std::advance(it, i);
        
        for (size_t j = 0; j < nodes.size(); ++j)
        {
            if ( nodes[j] == *it )
            {
                affected_nodes.erase(*it);
                --i;
                break;
            }
            
        }
        
    }
    
}


AbstractMove::AbstractMove( const AbstractMove &move ) : Move( move ),
    nodes( move.nodes ),
    affected_nodes( move.affected_nodes ),
    weight( move.weight ),
    auto_tuning( move.auto_tuning  ),
    num_tried_current_period( move.num_tried_current_period ),
    num_tried_total( move.num_tried_total )
{
    
    
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = nodes[i];
        
        // add myself to the set of moves
        the_node->addMove( this );
        
        // increase the DAG node reference count because we also have a pointer to it
        the_node->incrementReferenceCount();
        
    }
    
}



/**
 * Basic destructor doing nothing.
 */
AbstractMove::~AbstractMove( void )
{
    
    // clean up my pointers to the nodes
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = nodes[i];
        
        // add myself to the set of moves
        the_node->removeMove( this );
        
        // decrease the DAG node reference count because we also have a pointer to it
        if ( the_node->decrementReferenceCount() == 0 )
        {
            delete the_node;
        }
        
    }
    
}

/**
 * Overloaded assignment operator.
 * We need a deep copy of the operator.
 */
AbstractMove& AbstractMove::operator=(const RevBayesCore::AbstractMove &move)
{
    
    if ( this != &move )
    {
        // delegate
        Move::operator=(move);
        
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            // get the pointer to the current node
            DagNode* the_node = nodes[i];
            
            // add myself to the set of moves
            the_node->removeMove( this );
            
            // decrease the DAG node reference count because we also have a pointer to it
            if ( the_node->decrementReferenceCount() == 0 )
            {
                delete the_node;
            }
            
        }
        
        affected_nodes              = move.affected_nodes;
        nodes                       = move.nodes;
        num_tried_current_period    = move.num_tried_current_period;
        num_tried_total             = move.num_tried_total;
        
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            // get the pointer to the current node
            DagNode* the_node = nodes[i];
            
            // add myself to the set of moves
            the_node->addMove( this );
            
            // increase the DAG node reference count because we also have a pointer to it
            the_node->incrementReferenceCount();
            
        }
        
    }
    
    return *this;
}


/**
 * Add this node to our set of nodes.
 */
void AbstractMove::addNode( DagNode *n )
{
    // safety checking that the node is not NULL
    if ( n == NULL ) return;
    
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
    if ( exists == false )
    {
        nodes.push_back( n );
        
        // add myself to the set of moves
        n->addMove( this );
        
        // increment reference count
        n->incrementReferenceCount();
    }
    
}


/**
 * Automatically tune the move.
 * This will delegate the call to the derived class but only if autotuning was enabled.
 */
void AbstractMove::autoTune( void )
{
    
    // only call tuning if it is enabled for this move.
    if ( auto_tuning )
    {
        tune();
        
        // we need to reset the counters so that the next tuning only uses the new acceptance rate.
        resetCounters();
    }
    
}



/**
 * Decrement the counter for the number of tried attempts.
 */
void AbstractMove::decrementTriedCounter( void )
{
    // decrement the tries counter
    --num_tried_current_period;
    --num_tried_total;
    
}


/**
 * Get the set of affected nodes by a this move excluding the "core" nodes.
 *
 * \return The set of affected nodes.
 */
const RbOrderedSet<DagNode*>& AbstractMove::getAffectedNodes( void ) const
{
    
    return affected_nodes;
}


/**
 * Get the set of nodes on which this move is working on.
 *
 * \return The set of nodes.
 */
const std::vector<DagNode*>& AbstractMove::getDagNodes( void ) const
{
    
    return nodes;
}


/**
 * Get the number of how often the move has been used.
 *
 * \return    The number of tries.
 */
size_t AbstractMove::getNumberTriedCurrentPeriod( void ) const
{
    return num_tried_current_period;
}


/**
 * Get the number of how often the move has been used.
 *
 * \return    The number of tries.
 */
size_t AbstractMove::getNumberTriedTotal( void ) const
{
    return num_tried_total;
}


/**
 * Get the number of how often the move has been used.
 *
 * \return    The number of accepted moves.
 */
size_t AbstractMove::getNumberAcceptedCurrentPeriod( void ) const
{
    return num_tried_current_period;
}


/**
 * Get the number of how often the move has been used.
 *
 * \return    The number of accepted moves.
 */
size_t AbstractMove::getNumberAcceptedTotal( void ) const
{
    return num_tried_total;
}


/**
 * Get the update weight of how often the move should be used.
 *
 * \return    The update weight.
 */
double AbstractMove::getUpdateWeight( void ) const 
{
    return weight;
}


/**
 * Is the move active at the given generation?
 * The move itself will determine if it was switched off before or after
 * some given iteration.
 *
 * @param   generation    Generation number at which to check active status.
 */
bool AbstractMove::isActive(std::uint64_t generation) const
{
    return true;
}


/**
 * Dummy implementation of the hill-climbing move.
 */
void AbstractMove::performHillClimbingMove(double lHeat, double pHeat)
{
    
    throw RbException("A '" + getMoveName() + "' cannot be used for the hill-climbing algorithm." );
}



/**
 * Perform the move.
 * Here we store some info and delegate to performMove.
 */
void AbstractMove::performHillClimbingStep( double lHeat, double pHeat )
{
    // increment the tries counter
    ++num_tried_current_period;
    ++num_tried_total;
    
    // delegate to derived class
    performHillClimbingMove(lHeat, pHeat);
    
}



/**
 * Perform the move. 
 * Here we store some info and delegate to performMove.
 */
void AbstractMove::performMcmcStep( double prHeat, double lHeat, double pHeat )
{
    // increment the tries counter
    ++num_tried_current_period;
    ++num_tried_total;
    
    // delegate to derived class
    performMcmcMove(prHeat, lHeat, pHeat);
    
}




/**
 * Remove this node from our list of nodes.
 *
 * @param   node    The node to remove.
 */
void AbstractMove::removeNode( RevBayesCore::DagNode *node )
{
    
    for (std::vector<DagNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        if ( *it == node )
        {
            nodes.erase( it );
            break;
        }
    }
    
    // remove myself from this node
    node->removeMove( this );
    
    if ( node->decrementReferenceCount() == 0 )
    {
        delete node;
    }
    
}



/**
 * Reset the counters. Here we reset the counters stored in the base class
 * and delegate the call to the potentially overwritten method from the 
 * derived classes.
 *
 */
void AbstractMove::resetCounters( void )
{
    num_tried_current_period = 0;
    
    // delegate call
    resetMoveCounters();
}

/**
 * Reset the move counters. This method is only a hook for derived classes 
 * and here we provide only a dummy implementation.
 *
 */
void AbstractMove::resetMoveCounters( void )
{
    // dummy implementation
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     old_node     The old variable that needs to be replaced.
 * \param[in]     new_node     The new variable.
 */
void AbstractMove::swapNode(DagNode *old_node, DagNode *new_node)
{
    
    // find the old node
    
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // get the pointer to the current node
        DagNode* the_node = nodes[i];
        if ( the_node == old_node )
        {
            nodes[i] = new_node;
        }
        
    }
    
    // remove myself from the old node and add myself to the new node
    old_node->removeMove( this );
    new_node->addMove( this );
    
    // increment and decrement the reference counts
    new_node->incrementReferenceCount();
    if ( old_node->decrementReferenceCount() == 0 )
    {
        throw RbException("Memory leak in Metropolis-Hastings move. Please report this bug to Sebastian.");
    }
    
    affected_nodes.clear();
    
    for (std::vector<DagNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        // get the pointer to the current node
        DagNode* the_node = *it;
    
        // get the affected nodes if we would update this node
        // then we don't need to get the affected nodes every time again
        the_node->initiateGetAffectedNodes( affected_nodes );
        
    }
    
    // remove all "core" nodes from affectedNodes so their probabilities are not double-counted
    for (size_t i = 0; i < affected_nodes.size(); ++i)
    {
        RbOrderedSet<DagNode*>::iterator it = affected_nodes.begin();
        std::advance(it, i);
        
        for (size_t j = 0; j < nodes.size(); ++j)
        {
            if ( nodes[j] == *it )
            {
                affected_nodes.erase(*it);
                --i;
                break;
            }
            
        }
        
    }
    
    swapNodeInternal(old_node, new_node);
    
}


void AbstractMove::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_tried_current_period = na;
}


void AbstractMove::setNumberAcceptedTotal( size_t na )
{
    num_tried_total = na;
}


void AbstractMove::setNumberTriedCurrentPeriod( size_t nt )
{
    num_tried_current_period = nt;
}


void AbstractMove::setNumberTriedTotal( size_t nt )
{
    num_tried_total = nt;
}







