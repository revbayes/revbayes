#include "DagNode.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>

#include "Monitor.h"
#include "Move.h"
#include "RbException.h"
#include "RbOrderedSet.h"

using namespace RevBayesCore;

/**
 * Construct an empty DAG node, potentially
 * with a name (default is "").
 */
DagNode::DagNode( const std::string &n ) : Parallelizable(),
    children(),
    elementVar( false ),
    hidden( false ),
    monitors(),
    moves(),
    name( n ),
    touched_elements(),
    ref_count( 0 ),
    visit_flags( std::vector<bool>(5, false) )
{

}


/**
 * Copy constructor. We delegate the handling
 * of parents to derived classes with parents,
 * and they need to take care of the management
 * of those parents, removing them as a child
 * from those parents and deleting the parents
 * if necessary.
 */
DagNode::DagNode( const DagNode &n ) : Parallelizable( n ),
    children(),
    elementVar( n.elementVar ),
    hidden( n.hidden ),
    monitors( ),
    moves( ),
    name( n.name ),
    touched_elements( n.touched_elements ),
    ref_count( 0 ),
    visit_flags( n.visit_flags )
{

}


/**
 * Destructor. The destructor should never be called when the reference count
 * is larger than 0, or when we have children. All children should increase the
 * reference count, so there should be no children left when the reference count
 * is 0. Other Rev objects may also increase the reference count to protect a
 * DAG node from being deleted before they die, so the reference count can be
 * larger than the number of children, but never smaller.
 *
 * Deletion of parents is delegated to the destructors of derived classes with
 * parents.
 */
DagNode::~DagNode( void )
{

}


/**
 * Assignment operator. Note that parent management is delegated to
 * derived classes that actually have parents, so their assignment
 * operators need to deal with parent replacement at assignment.
 *
 * Note that children do not need to replace or remove this node
 * as their parent, because we remain alive and their parent
 * pointer will continue to be relevant.
 */
DagNode& DagNode::operator=(const DagNode &d)
{
    Parallelizable::operator=( d );

    if ( &d != this )
    {
        name             = d.name;
        elementVar       = d.elementVar;
        hidden           = d.hidden;
        touched_elements = d.touched_elements;
        visit_flags      = d.visit_flags;
    }

    return *this;
}


/**
 * Add a new child node to this node.
 * Since we store the children in a set we don't need to worry about duplicates.
 *
 * Note, the caller also needs to increment the reference count to this node.
 */
void DagNode::addChild(DagNode *child) const
{

    // only if the child is not NULL and isn't in our vector yet
    if ( child != NULL )
    {
        std::vector<DagNode*>::const_iterator pos = std::find(children.begin(), children.end(), child);
        if ( pos == children.end() )
        {
            children.push_back( child );
        }

    }

}


/**
 * Add a new monitor which monitors this node.
 * We only keep these pointers to notify the monitor if we are replaced.
 *
 *
 * Note, the caller also needs to increment the reference count to this node.
 */
void DagNode::addMonitor(Monitor *m)
{
    // only if the child is not NULL and isn't in our vector yet
    if ( m != NULL )
    {
        std::vector<Monitor *>::const_iterator pos = std::find(monitors.begin(), monitors.end(), m);
        if ( pos == monitors.end() )
        {
            monitors.push_back( m );
        }

    }

}


/**
 * Add a new move which update this node.
 * We only keep these pointers to notify the monitor if we are replaced.
 *
 *
 * Note, the caller also needs to increment the reference count to this node.
 */
void DagNode::addMove(Move *m)
{
    // only if the child is not NULL and isn't in our vector yet
    if ( m != NULL )
    {
        std::vector<Move *>::const_iterator pos = std::find(moves.begin(), moves.end(), m);
        if ( pos == moves.end() )
        {
            moves.push_back( m );
        }

    }

}


/* Add the index of an element that has been touched */
void DagNode::addTouchedElementIndex(size_t i)
{

    touched_elements.insert( i );

}


void DagNode::clearVisitFlag( const size_t& flagType )
{

    RbOrderedSet<DagNode*> descendants;
    findUniqueDescendantsWithFlag(descendants, flagType);

    // Clear the designated flagType from all descedants (including node calling this)
    // Also clear the flags we just flagged to keep descedant searching fast
    for (RbOrderedSet<DagNode*>::const_iterator it=descendants.begin(); it!=descendants.end(); ++it)
    {
        DagNode *the_descendant = *it;
        the_descendant->visit_flags[FIND_FLAG] = false;
        the_descendant->visit_flags[flagType] = false;
    }


    return;
}


void DagNode::clearTouchedElementIndices( void )
{

    touched_elements.clear();

}


/** Clone the graph downstream from this node: clone children */
DagNode* DagNode::cloneDownstreamDag( std::map<const DagNode*, DagNode* >& newNodes ) const
{

    if ( newNodes.find( this ) != newNodes.end() )
    {
        return ( newNodes[ this ] );
    }
    // Get pristine copy
    DagNode* copy = clone();

    // Add this node to the map
    newNodes[ this ] = copy;

    /* Make sure the children clone themselves */
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = (*it)->cloneDownstreamDag( newNodes );
        child->swapParent(this, copy);
    }

    return copy;
}


/**
 * Decrement the reference count and return it.
 */
size_t DagNode::decrementReferenceCount( void ) const
{

    --ref_count;

    return ref_count;
}


void DagNode::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const
{

    if ( n == "lnProbability" )
    {
        rv = const_cast<DagNode *>(this)->getLnProbability();
    }
    else if ( n == "probability" )
    {
        rv = std::exp(const_cast<DagNode *>(this)->getLnProbability());
    }
    else
    {
        throw RbException() << "A DAG node does not have a member method called '" << n << "'.";
    }

}

/*
 * finds all descendants without redundant node visitation
 */
void DagNode::findUniqueDescendants(RbOrderedSet<DagNode *>& descendants)
{
    visit_flags[FIND_FLAG] = true;

    // add self to descendant list
    descendants.insert(this);

    // recurse across node's children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        // if child is not in descedant list, recurse from child's position
        // if ( descendants.find( *it ) == descendants.end() )
        if (child->visit_flags[FIND_FLAG] == false )
        {
            child->findUniqueDescendants( descendants );
        }
    }
}


/*
 * finds all descendants without redundant node visitation
 */
void DagNode::findUniqueDescendantsWithFlag(RbOrderedSet<DagNode *>& descendants, const size_t flagType)
{
    visit_flags[FIND_FLAG] = true;

    // add self to descendant list
    descendants.insert(this);

    // recurse across node's children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        // if child is not in descedant list, recurse from child's position
        // if ( descendants.find( *it ) == descendants.end() )
        if ( child->visit_flags[FIND_FLAG] == false && child->visit_flags[flagType] == true )
        {
            child->findUniqueDescendantsWithFlag( descendants, flagType );
        }
    }
}


/**
 * Get all affected nodes this DAGNode.
 * This means we call getAffected() of all children. getAffected() is pure virtual.
 */
void DagNode::getAffectedNodes(RbOrderedSet<DagNode *> &affected) const
{
    visit_flags[AFFECTED_FLAG] = true;
    // get all my affected children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        if (child->visit_flags[AFFECTED_FLAG] == false)
        {
            child->getAffected(affected, this);
        }
    }
}


/**
 * Get a const reference to the local set of children of this node.
 */
const std::vector<DagNode*>& DagNode::getChildren( void ) const
{
    return children;
}


/**
 * Get the type of the DAG node as a string.
 */
DagNode::DagNodeTypes DagNode::getDagNodeType( void ) const
{

    return type;
}


/**
 * Get the distribution of this node. Only stochastic nodes have a distribution and thus we throw an error here.
 */
const Distribution& DagNode::getDistribution( void ) const
{
    throw RbException("Only stochastic nodes have a distribution.");
}


/**
 * Get the distribution of this node. Only stochastic nodes have a distribution and thus we throw an error here.
 */
Distribution& DagNode::getDistribution( void )
{
    throw RbException("Only stochastic nodes have a distribution.");
}


/**
 * Get the first child of this node.
 * Here we simply return a pointer to the first element stored in the set of children.
 */
DagNode* DagNode::getFirstChild( void ) const
{

    return *children.begin();
}


std::vector<double> DagNode::getMixtureProbabilities(void) const
{
    return std::vector<double>(1,1.0);
}



/**
 * Get a const reference to our local monitor set.
 */
const std::vector<Monitor*>& DagNode::getMonitors( void ) const
{
    // return the local object
    return monitors;
}


/**
 * Get a const reference to our local moves set.
 */
const std::vector<Move*>& DagNode::getMoves( void ) const
{
    // return the local object
    return moves;
}


/**
 * Get a const reference to the variable name.
 */
const std::string& DagNode::getName( void ) const
{

    return name;
}


/**
 * Get the number of children for this DAG node.
 */
size_t DagNode::getNumberOfChildren( void ) const
{

    return children.size();
}


size_t DagNode::getNumberOfMixtureElements(void) const
{
    return 0;
}


/**
 * To allow calls to getParents even on classes without parents,
 * we simply return an empty set by default.
 */
std::vector<const DagNode*> DagNode::getParents( void ) const
{

    return std::vector<const DagNode*>();
}


/**
 * Get the printable children by filling the set of DAG nodes with the printable children.
 * This method will skip hidden variables and instead replace the hidden variables by their children,
 * thus calling this function recursively for the hidden children.
 * Hidden children may be converter nodes or similar.
 */
void DagNode::getPrintableChildren(std::vector<DagNode *> &c) const
{

    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        if ( child->isHidden() == false )
        {
            // just insert this child (if we haven't already)
            std::vector<DagNode *>::const_iterator pos = std::find(c.begin(), c.end(), child);
            if ( pos == c.end() )
            {
                c.push_back( child );
            }

        }
        else
        {
            // do not add this child but all the children below because we omit this node
            child->getPrintableChildren( c );
        }

    }

}



/**
 * Get the printable parents by filling the set of DAG nodes with the printable parents.
 * This method will skip hidden variables and instead replace the hidden variables by their parents,
 * thus calling this function recursively for the hidden parents.
 * Hidden parents may be converter nodes or similar.
 */
void DagNode::getPrintableParents(std::vector<const DagNode *> &p) const
{

    std::vector<const DagNode*> parents = getParents();
    for (std::vector<const DagNode*>::const_iterator it=parents.begin(); it!=parents.end(); ++it)
    {
        const DagNode *parent = *it;
        if ( parent->isHidden() == false )
        {
            // just insert this child (if we haven't already)
            std::vector<const DagNode *>::const_iterator pos = std::find(p.begin(), p.end(), parent);
            if ( pos == p.end() )
            {
                p.push_back( parent );
            }

        }
        else
        {
            // do not add this child but all the children below because we omit this node
            parent->getPrintableParents( p );
        }

    }

}


/**
 * Get the reference count.
 */
size_t DagNode::getReferenceCount( void ) const
{

    return ref_count;
}


/* Get the indices of all touched elements */
const std::set<size_t>& DagNode::getTouchedElementIndices( void ) const
{

    return touched_elements;
}

bool DagNode::getVisitFlag( const size_t flagType ) const
{
    return visit_flags[flagType];
}


/**
 * Increment the reference count.
 */
void DagNode::incrementReferenceCount( void ) const
{

    ++ref_count;

}


/**
 * Begins a getAffectedNodes() recursion then clears visited flags
 */
void DagNode::initiateGetAffectedNodes(RbOrderedSet<DagNode *> &affected)
{

    // begin recursion
    getAffectedNodes( affected );

    // clear visit flags
    const size_t& flag_type = AFFECTED_FLAG;
    clearVisitFlag(flag_type);
}


/**
 * This function returns true if the DAG node is
 * a node that is modifiable by the user through
 * assignment. This occurs when the node is named
 * or at least one of the parents is assignable.
 *
 * @todo The current code really belongs to the
 *       language layer but it is convenient to
 *       implement it here. See if it is possible
 *       to move it to the language layer.
 */
bool DagNode::isAssignable( void ) const
{
    const std::vector<const DagNode*>& parents = getParents();

    if ( getName() != "" )
    {
        return true;
    }

    for (std::vector<const DagNode*>::const_iterator it=parents.begin(); it!=parents.end(); ++it)
    {
        const DagNode *parent = *it;
        if ( parent->isAssignable() )
        {
            return true;
        }
    }

    return false;
}


bool DagNode::isClamped( void ) const
{

    return false;
}


bool DagNode::isConstant( void ) const
{

    return false;
}

bool DagNode::isElementVariable( void ) const
{

    return elementVar;
}

bool DagNode::isHidden( void ) const
{

    return hidden;
}


bool DagNode::isIntegratedOut( void ) const
{
    return false;
}


/**
 * Is this variable a simple numeric variable?
 * This is asked for example by the model monitor that only wants to monitor simple numeric variable
 * because all others (e.g. trees and vectors/matrices) cannot be read by Tracer.
 */
bool DagNode::isSimpleNumeric( void ) const
{
    return false;
}


bool DagNode::isStochastic( void ) const
{
    return false;
}

/**
 * Keep the value of the node.
 * This function delegates the call to keepMe() and calls keepAffected() too.
 */
void DagNode::keep(void)
{

    // keep myself first
    keepMe( this );

    // next, keep all my children
    keepAffected();

    // clear visit flags
    const size_t &flag_type = KEEP_FLAG;
    clearVisitFlag(flag_type);

}


/**
 * Tell affected variable nodes to keep the current value.
 */
void DagNode::keepAffected()
{
    visit_flags[KEEP_FLAG] = true;

    // keep all my children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        if (child->visit_flags[KEEP_FLAG] == false)
        {
            child->keepMe( this );
        }
    }

}


/** Print children. We assume indent is from line 2 (hanging indent). Line length is lineLen. */
void DagNode::printChildren( std::ostream& o, size_t indent, size_t lineLen, bool verbose ) const
{
    std::string pad;
    for ( size_t i = 0; i < indent; ++i )
    {
        pad.push_back( ' ' );
    }

    o << "[ ";
    pad += "  ";
    indent += 2;

    size_t currentLength = indent + 2;
    std::ostringstream s;

    // create my own copy of the pointers to the children
    std::vector<DagNode*> printableChildren = children;

    // replace the children that should not be printed
    if ( verbose == false )
    {

        printableChildren.clear();
        getPrintableChildren( printableChildren );

    }

    std::vector<DagNode*>::const_iterator it;
    size_t i = 0;
    for ( i = 0, it = printableChildren.begin(); it != printableChildren.end(); ++it, ++i )
    {
        std::ostringstream s;
        if ( (*it)->getName() == "" )
        {
            s << "<" << (*it) << ">";
        }
        else
        {
            if ( verbose == true )
            {
            s << (*it)->getName() << " <" << (*it) << ">";
            }
            else
            {
                s << (*it)->getName();
            }
        }

        if ( printableChildren.size() - i > 1 )
        {
            s << ", ";
        }

        if ( s.str().size() + currentLength > lineLen )
        {
            o << std::endl << pad;
            currentLength = pad.size();
        }
        o << s.str();
        currentLength += s.str().size();
        s.str("");
    }

    o << " ]";
}


/**
 * Print parents. This function does not access the parents directly, but
 * through the getParents function. This is to allow derived classes to
 * handle parents in a different way than done in this class.
 *
 * Note the use of a const reference here to potentially avoid one copy
 * operation on the set of DAG nodes.
 */
void DagNode::printParents( std::ostream& o, size_t indent, size_t lineLen, bool verbose ) const
{

    std::string pad;
    for ( size_t i = 0; i < indent; ++i )
    {
        pad.push_back( ' ' );
    }

    o << "[ ";
    pad += "  ";
    indent += 2;

    size_t currentLength = indent + 2;
    std::ostringstream s;

    // create my own copy of the pointers to the children
    std::vector<const DagNode*> printableParents = getParents();

    // replace the children that should not be printed
    if ( verbose == false )
    {

        printableParents.clear();
        getPrintableParents( printableParents );

    }

    std::vector<const DagNode*>::const_iterator it;
    size_t i = 0;

    for ( i = 0, it = printableParents.begin(); it != printableParents.end(); ++it, ++i )
    {
        if ( (*it)->getName() == "" )
        {
            s << "<" << (*it) << ">";
        }
        else
        {
            if ( verbose == true )
            {
                s << (*it)->getName() << " <" << (*it) << ">";
            }
            else
            {
                s << (*it)->getName();
            }
        }

        if ( printableParents.size() - i > 1 )
        {
            s << ", ";
        }

        if ( s.str().size() + currentLength > lineLen )
        {
            o << std::endl << pad;
            currentLength = pad.size();
        }
        o << s.str();
        currentLength += s.str().size();
        s.str("");
    }
    o << " ]";
}


/**
 * By default we do not need to do anything when re-initializiating.
 */
void DagNode::reInitialized( void )
{

    reInitializeMe();
    reInitializeAffected();

    // clear visit flags
    const size_t &flag_type = REINITIALIZE_FLAG;
    clearVisitFlag(flag_type);

}


/**
 * By default we do not need to do anything when re-initializiating.
 */
void DagNode::reInitializeAffected( void )
{

    visit_flags[REINITIALIZE_FLAG] = true;

    // next, reInitialize all my children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        if (child->visit_flags[REINITIALIZE_FLAG] == false)
        {
            child->reInitializeMe();
        }
    }

}


/**
 * By default we do not need to do anything when re-initializiating.
 */
void DagNode::reInitializeMe( void )
{
    // nothing to do
}


/**
 * Remove the DAG node from our set of children.
 *
 * Note, the caller also needs to decrement the reference count to this node.
 */
void DagNode::removeChild(DagNode *child) const
{

    // test if we even have this node as a child
    std::vector<DagNode *>::const_iterator it = std::find( children.begin(), children.end(), child );
    if ( it != children.end() )
    {
        children.erase( it );

        // we do not own our children! See addChildNode for explanation

    }

}





/**
 * Remove the monitor from our set of monitors.
 *
 * Note, the caller also needs to decrement the reference count to this node.
 */
void DagNode::removeMonitor(Monitor *m)
{

    // test if we even have this monitor
    std::vector<Monitor*>::const_iterator it = std::find( monitors.begin(), monitors.end(), m );
    if ( it != monitors.end() )
    {
        monitors.erase( it );

        // we do not own our monitors! See addMonitor for explanation

    }

}


/**
 * Remove the move from our set of moves.
 *
 * Note, the caller also needs to decrement the reference count to this node.
 */
void DagNode::removeMove(Move *m)
{

    // test if we even have this move
    std::vector<Move*>::const_iterator it = std::find( moves.begin(), moves.end(), m );
    if ( it != moves.end() )
    {
        moves.erase( it );

        // we do not own our moves! See addMoves for explanation

    }

}


/**
 * Replace the DAG node.
 * We call swap parent for all children so that they get a new parent. We do not change the parents.
 */
void DagNode::replace( DagNode *n )
{

    // replace myself for all my monitor
    while ( monitors.empty() == false )
    {
        Monitor *theMonitor = *monitors.begin();
        theMonitor->swapNode( this, n);
    }

    // replace myself for all my moves
    while ( moves.empty() == false )
    {
        Move *the_move = *moves.begin();
        the_move->swapNode( this, n);
    }

    // replace myself at all my children
    while ( children.empty() == false )
    {
        DagNode *child = *children.begin();
        child->swapParent( this, n);
    }

}


/**
 * Restore this DAGNode.
 * This means we call restoreMe() and restoreAffected(). There is a default implementation of restoreAffected()
 * which call restoreMe() of all children of this node. restoreMe() is pure virtual.
 */
void DagNode::restore(void)
{

    // first restore myself
    restoreMe( this );

    // next, restore all my children
    restoreAffected();

    // clear visit flags
    const size_t &flag_type = RESTORE_FLAG;
    clearVisitFlag(flag_type);

}


/**
 * Restore all affected nodes this DAGNode.
 * This means we call restoreMe() of all children. restoreMe() is pure virtual.
 */
void DagNode::restoreAffected(void)
{

    visit_flags[RESTORE_FLAG] = true;

    // keep all my children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        if (child->visit_flags[RESTORE_FLAG] == false)
        {
            child->restoreMe( this );
        }
    }
}


void DagNode::setElementVariable(bool tf)
{
    elementVar = tf;
}


void DagNode::setHidden(bool tf)
{
    hidden = tf;
}


void DagNode::setIntegrationIndex( size_t i )
{
    // nothing to do in the base class
}



void DagNode::setName(std::string const &n)
{
    // set the internal value
    name = n;

}


void DagNode::setParentNamePrefix(const std::string &p)
{
    std::vector<const DagNode*> parents = getParents();
    for (std::vector<const DagNode*>::const_iterator it=parents.begin(); it!=parents.end(); ++it)
    {
        const DagNode *parent = *it;
        if ( parent->getName().size() > 0 && parent->getName()[0] == '.' )
        {
            DagNode *parentNode = const_cast<DagNode *>( parent );
            // just insert this child
            parentNode->setName( p + parentNode->getName() );
            parentNode->setParentNamePrefix( parentNode->getName() );
        }
    }

}


void DagNode::setPriorOnly(bool tf)
{
    throw RbException()<<"setPriorOnly("<<tf<<"): node is not stochastic!";

}

void DagNode::setVisitFlag(bool tf, const size_t flagType)
{
    visit_flags[flagType] = tf;
}

/**
 * Swap parent node. We delegate this task to derived DAG node classes
 * with parents. Here we just throw an error in case a derived class
 * with parents forgets to override this function.
 */
void DagNode::swapParent( const DagNode *oldParent, const DagNode *newParent )
{

    throw RbException( "This DAG node does not have any parents" );

}


/**
 * Touch the DAG node.
 *
 * This function should be called if the value of the variable has changed or if you want this node to be reevaluated.
 * The function will automatically call the touchMe() which is implemented differently in the different DAG node types.
 *
 * Since the DAG node was touched and possibly changed, we tell affected DAG nodes that they too have been touched
 * and need to update their value.
 */
void DagNode::touch(bool touchAll)
{
    // first touch myself
    touchMe( this, touchAll );

    // next, touch all my children
    touchAffected( touchAll );
}


/**
 * Tell affected variable nodes to touch themselves (i.e. that they've been touched).
 */
void DagNode::touchAffected(bool touchAll)
{
    // touch all my children
    for (std::vector<DagNode*>::const_iterator it=children.begin(); it!=children.end(); ++it)
    {
        DagNode *child = *it;
        child->touchMe( this, touchAll );
    }
}

double DagNode::getPrevLnProbability(void) const
{
    throw RbException()<<"getPrevLnProbability: not a stochastic node!";
}
