#include "TraceTree.h"

using namespace RevBayesCore;


/*
 * TraceTree constructor
 */
TraceTree::TraceTree( bool c ) : TreeSummary(this, c)
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
