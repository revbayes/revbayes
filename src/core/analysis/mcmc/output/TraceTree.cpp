#include "TraceTree.h"

using namespace RevBayesCore;


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
    TreeSummary(this, c),
    clock( c ),
    rooted( c )
{
}


/*
 * TraceTree copy constructor
 */
TraceTree::TraceTree(const TraceTree& t ) : TreeSummary(this, t.isClock())
{
    *this = t;

    traces.clear();
    traces.push_back(this);
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
    return (TreeSummary::isCoveredInInterval(v,size,verbose) ? 0 : -1);
}

int TraceTree::isCoveredInInterval(const Tree &t, double size, bool verbose)
{
    return (TreeSummary::isCoveredInInterval(t,size,verbose) ? 0 : -1);
}

