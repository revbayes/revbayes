#include "TraceTree.h"
#include "TreeSummary.h"

using namespace RevBayesCore;


const TreeSummary& TraceTree::summary() const
{
    return *the_summary;
}

TreeSummary& TraceTree::summary()
{
    return *the_summary;
}

void TraceTree::doAddObject(Tree&& t)
{
    if (isRooted() and not t.isRooted())
        throw RbException()<<"Adding unrooted tree to rooted trace";
    else if (not isRooted() and t.isRooted())
        throw RbException()<<"Adding rooted tree to unrooted trace";

    return Trace<Tree>::doAddObject( std::move(t) );
}

/*
 * TraceTree constructor
 */
TraceTree::TraceTree( bool c ) :
    clock( c ),
    rooted( c ),
    the_summary( std::make_unique<TreeSummary>(this, clock) )
{
}


/*
 * TraceTree copy constructor and operator=.
 * Hopefully we can eventually delete these.
 */

TraceTree::TraceTree(const TraceTree& t)
{
    operator=(t);
}

TraceTree&  TraceTree::operator=(const TraceTree& t)
{
    Trace<Tree>::operator=(t);

    clock = t.clock;
    rooted = t.rooted;

    the_summary = std::make_unique<TreeSummary>(this, clock);

    if (t.summary().hasOutgroup())
        summary().setOutgroup(t.summary().getOutgroup());

    return *this;
}

bool TraceTree::isRooted() const
{
    return rooted;
}

bool TraceTree::isClock() const
{
    return clock;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */

TraceTree* TraceTree::clone(void) const
{
    
    return new TraceTree(*this);
}

int TraceTree::isCoveredInInterval(const std::string &v, double size, bool verbose)
{
    return (summary().isCoveredInInterval(v,size,verbose) ? 0 : -1);
}

int TraceTree::isCoveredInInterval(const Tree &t, double size, bool verbose)
{
    return (summary().isCoveredInInterval(t,size,verbose) ? 0 : -1);
}

