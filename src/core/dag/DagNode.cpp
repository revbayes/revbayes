#include "DagNode.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Monitor.h"
#include "Move.h"
#include "RbException.h"
#include "RbOrderedSet.h"

using namespace RevBayesCore;


/**
 * Online cycle detection via Pearce-Kelly topological-order maintenance: see
 *
 * Pearce DJ, Kelly PHJ. 2007. A dynamic topological sort algorithm for directed acyclic graphs. Journal of Experimental
 * Algorithmics (JEA). 11:1.7--es. doi:10.1145/1187436.121059.
 *
 * Each `DagNode` stores an ordinal value `topo_ord` (its rank in the current topological order). The defining invariant is:
 *
 *     for every edge (parent -> child) in the DAG,
 *         parent->topo_ord < child->topo_ord.
 *
 * Ordinals are assigned monotonically by the static counter `next_topo_ord`: each constructor (default and copy) grabs the
 * current value and increments. Since every existing node has a strictly smaller ordinal than the value a brand-new node
 * receives, the invariant automatically holds for every edge of the form "existing parent -> brand-new child". The only
 * operation that can violate the invariant is reattaching a previously existing node above another previously existing node.
 * When `addChild()` notices `topo_ord >= child->topo_ord`, it delegates to `pkAddEdge()`, which:
 *
 *   1. Runs a forward depth-first search from `child`, restricted to nodes whose ordinal is <= ord(parent). This is the
 *      "forward delta", deltaF, of nodes that need to be promoted so they end up above `parent` in the new order. Reaching
 *      `parent` during this search would mean the new edge closes a cycle, so we throw without mutating any ordinal.
 *   2. Runs a backward depth-first search from `parent`, restricted to nodes whose ordinal is >= ord(child). This is the
 *      "backward delta", deltaB, of nodes that need to be demoted so they end up below `child`.
 *   3. Reassigns the original ordinals of the union of deltaB and deltaF so that every node in deltaB precedes every node
 *      in deltaF in the new ordering. The original ordinals are recycled in place, so the ordinal space stays compact within
 *      the affected band. Each delta is collected in ascending ordinal order via a post-order depth-first walk with neighbors
 *      sorted at each step, and the multiset of ordinals is merged using std::merge in O(|deltaB| + |deltaF|).
 *
 * `dependsOn(node, target)` then becomes cheap: an ancestor must have a strictly smaller ordinal than its descendant,
 * so `ord(target) >= ord(node)` answers "no" in O(1). Otherwise we walk up the parents of `node`, but pruned to the band
 * `ord >= ord(target)`, which keeps the search inside a small region of the DAG instead of the whole upward cone.
 *
 * Edge removal by `removeChild()` preserves the invariant trivially (it can only widen the set of valid topological orderings),
 * so it needs no maintenance.
 */

// Graph mutation in RevBayes is single-threaded; the global counter is therefore a plain size_t (no need for atomics).
std::size_t DagNode::next_topo_ord = 0;


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
    topo_ord( next_topo_ord++ ),
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
    topo_ord( next_topo_ord++ ),
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
 * Add a new child node to this node. Since we store the children in a set, we don't need to worry about duplicates.
 *
 * If the new edge violates the Pearce-Kelly topological-order invariant (when `this` was created after `child` and is now being
 * attached above it, e.g., through `swapParent()`), we delegate to `pkAddEdge` to repair it. If `pkAddEdge` detects that the
 * new edge would close a cycle, it throws; we then roll the children vector back so that the DAG is left untouched.
 *
 * Note, the caller also needs to increment the reference count to this node.
 */
void DagNode::addChild(DagNode *child) const
{

    // only if the child is not NULL and isn't in our vector yet
    if ( child == NULL )
    {
        return;
    }

    std::vector<DagNode*>::const_iterator pos = std::find(children.begin(), children.end(), child);
    if ( pos != children.end() )
    {
        return;
    }

    children.push_back( child );

    // Repair the topological ordering if the new edge violates the invariant.
    if ( topo_ord >= child->topo_ord )
    {
        try
        {
            pkAddEdge( this, child );
        }
        catch ( ... )
        {
            // Roll back the edge so the graph state matches what the caller saw
            // before the failed insertion.
            children.pop_back();
            throw;
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


/**
 * Find out whether `node`, or any of its ancestors, equals `target`.
 *
 * The DAG keeps a topological order via the Pearce-Kelly machinery, so if `target` is an ancestor of `node`, then
 * ord(target) < ord(node). The contrapositive lets us reject no-cycle case in O(1): if ord(target) >= ord(node),
 * `target` cannot be an ancestor.
 *
 * In the slow path, we walk upward from `node` along parent edges, but only visit nodes whose ordinal lies in
 * [ord(target), ord(node)]. This keeps the search confined to the affected band of the DAG instead of traversing
 * the entire upward cone, which would exhibit quadratic total cost.
 */
 bool DagNode::dependsOn(const DagNode* node, const DagNode* target)
 {
     // Fast path: an ancestor must have a strictly smaller ordinal than its descendant, so anything at or above `node`
     // in the ordering cannot possibly be an ancestor of `node`.
     if ( target->topo_ord >= node->topo_ord )
     {
         return false;
     }
 
     const std::size_t lb = target->topo_ord;
 
     std::unordered_set<const DagNode*> visited;
     std::vector<const DagNode*> stack;
     stack.push_back( node );
 
     while ( not stack.empty() )
     {
         const DagNode* current = stack.back();
         stack.pop_back();
 
         if ( visited.count( current ) ) continue;
         visited.insert( current );
 
         if ( current == target ) return true;
 
         // Only follow parents that could conceivably still reach `target`. Any parent with ord < lb is too far up the order
         // to be `target` or to have `target` among its own ancestors.
         for ( const DagNode* parent : current->getParents() )
         {
             if ( parent->topo_ord >= lb && visited.count( parent ) == 0 )
             {
                 stack.push_back( parent );
             }
         }
     }
 
     return false;
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
 * Get a name we can display when throwing a warning that a cycle has been detected.
 */
std::string DagNode::getNodeDisplayName(const DagNode* n)
{
    if ( n == nullptr ) return "<null>";
    const std::string& nm = n->getName();
    if ( !nm.empty() ) return nm;
    switch ( n->getDagNodeType() )
    {
        case DagNode::CONSTANT:      return "<unnamed> (constant node)";
        case DagNode::DETERMINISTIC: return "<unnamed> (deterministic node)";
        case DagNode::STOCHASTIC:    return "<unnamed> (stochastic node)";
    }
    return "<unnamed>";
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


double DagNode::getPrevLnProbability(void) const
{
    throw RbException()<<"getPrevLnProbability: not a stochastic node!";
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


/**
 * Get the position of this node in the maintained topological ordering of the DAG. The Pearce-Kelly invariant guarantees that
 * for every edge `parent -> child`, `parent->getTopologicalOrder() < child->getTopologicalOrder()`.
 */
 std::size_t DagNode::getTopologicalOrder( void ) const
 {
     return topo_ord;
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

bool DagNode::isIgnoredData( void ) const
{

    return false;
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


/**
 * Pearce-Kelly online edge insertion: repair the topological-order invariant after the edge `parent -> child` has been added
 * with `ord(parent) >= ord(child)`, or throw if the new edge would close a cycle.
 *
 *     deltaF (forward)  = nodes reachable from `child` via children with ord <= ord(parent)
 *     deltaB (backward) = nodes reachable from `parent` via parents with ord >= ord(child)
 *
 * The forward search doubles as cycle detection: if it ever reaches `parent`, the new edge would close a cycle, and we throw
 * without mutating any ordinal.
 *
 * Each delta must be listed in ascending ordinal order before reassignment. We avoid sorting the whole delta by doing an
 * explicit post-order traversal and, at each node, sorting only that node's filtered neighbors by ordinal before recursing (so
 * the smallest-ordinal child branches finish first). The forward walk records finishes in descending ordinal order, so we
 * reverse once; the backward walk's finish order is already ascending. The pooled ordinals are then merged in linear time.
 */
void DagNode::pkAddEdge( const DagNode* parent, DagNode* child )
{
    const std::size_t lb = child->topo_ord;
    const std::size_t ub = parent->topo_ord;

    auto by_ord = []( DagNode* a, DagNode* b ) { return a->topo_ord < b->topo_ord; };

    // Forward post-order from child along children (ord <= ub). Entering `parent` <=> cycle.
    std::vector<DagNode*> deltaF;
    std::unordered_set<DagNode*> visitedF;
    {
        std::vector<std::pair<DagNode*, bool> > stack;
        stack.push_back( std::make_pair( child, false ) );
        std::vector<DagNode*> next_children;
        while ( not stack.empty() )
        {
            std::pair<DagNode*, bool> frame = stack.back();
            stack.pop_back();
            DagNode* z = frame.first;
            const bool is_exit = frame.second;

            if ( is_exit )
            {
                deltaF.push_back( z );
                continue;
            }

            if ( visitedF.count( z ) ) continue;

            if ( z == parent )
            {
                throw RbException( "Cannot add this edge to the DAG: it would create a cycle." );
            }

            visitedF.insert( z );
            stack.push_back( std::make_pair( z, true ) );

            next_children.clear();
            for ( DagNode* d : z->getChildren() )
            {
                if ( d->topo_ord <= ub )
                {
                    next_children.push_back( d );
                }
            }
            std::sort( next_children.begin(), next_children.end(), by_ord );
            for ( std::vector<DagNode*>::reverse_iterator it = next_children.rbegin(); it != next_children.rend(); ++it )
            {
                stack.push_back( std::make_pair( *it, false ) );
            }
        }
        std::reverse( deltaF.begin(), deltaF.end() );
    }

    // Backward post-order from parent along parents (ord >= lb). Finish order is ascending ordinal.
    std::vector<DagNode*> deltaB;
    std::unordered_set<DagNode*> visitedB;
    {
        std::vector<std::pair<DagNode*, bool> > stack;
        stack.push_back( std::make_pair( const_cast<DagNode*>( parent ), false ) );
        std::vector<DagNode*> next_parents;
        while ( not stack.empty() )
        {
            std::pair<DagNode*, bool> frame = stack.back();
            stack.pop_back();
            DagNode* z = frame.first;
            const bool is_exit = frame.second;

            if ( is_exit )
            {
                deltaB.push_back( z );
                continue;
            }

            if ( visitedB.count( z ) ) continue;

            visitedB.insert( z );
            stack.push_back( std::make_pair( z, true ) );

            next_parents.clear();
            for ( const DagNode* p : z->getParents() )
            {
                DagNode* pm = const_cast<DagNode*>( p );
                if ( pm->topo_ord >= lb )
                {
                    next_parents.push_back( pm );
                }
            }
            std::sort( next_parents.begin(), next_parents.end(), by_ord );
            for ( std::vector<DagNode*>::reverse_iterator it = next_parents.rbegin(); it != next_parents.rend(); ++it )
            {
                stack.push_back( std::make_pair( *it, false ) );
            }
        }
    }

    std::vector<std::size_t> ordsB;
    std::vector<std::size_t> ordsF;
    ordsB.reserve( deltaB.size() );
    ordsF.reserve( deltaF.size() );
    for ( DagNode* z : deltaB ) ordsB.push_back( z->topo_ord );
    for ( DagNode* z : deltaF ) ordsF.push_back( z->topo_ord );

    std::vector<std::size_t> pool( deltaB.size() + deltaF.size() );
    std::merge( ordsB.begin(), ordsB.end(), ordsF.begin(), ordsF.end(), pool.begin() );

    std::size_t i = 0;
    for ( DagNode* z : deltaB ) z->topo_ord = pool[i++];
    for ( DagNode* z : deltaF ) z->topo_ord = pool[i++];
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

    // Before modifying anything, verify the replacement would not create a cycle.
    // After replacement, all of this node's children become n's children.
    // If n already depends (transitively) on this node, that dependency chain
    // would close into a loop.
    if ( n != NULL && dependsOn(n, this) )
    {
        throw RbException( "Cannot complete this assignment: '" + getNodeDisplayName( n ) + "' depends on '" + getNodeDisplayName( this ) + "', so replacing one with the other would create a circular dependency." );
    }

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


void DagNode::setIgnoreData(bool tf)
{
    throw RbException()<<"Error: can't ignore data at node '"<<getName()<<"' because it is not a stochastic node!";

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
    for (DagNode* child: children)
    {
        child->touchMe( this, touchAll );
    }
}

