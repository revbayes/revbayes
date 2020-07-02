#ifndef TraceTree_H
#define TraceTree_H

#include "Clade.h"
#include "Trace.h"
#include "Tree.h"
#include "TreeSummary.h"

namespace RevBayesCore {

    class TraceTree : public Trace<Tree>, public TreeSummary {

        public:

        /*
         * Declaration of the TreeTrace class
         */
        TraceTree( bool c = true );
        TraceTree(const TraceTree& t );
        virtual ~TraceTree(){}

        TraceTree*                                 clone(void) const;

        int                                        isCoveredInInterval(const std::string &v, double size, bool verbose){ return (TreeSummary::isCoveredInInterval(v,size,verbose) ? 0 : -1); };
        int                                        isCoveredInInterval(const Tree &t, double size, bool verbose){ return (TreeSummary::isCoveredInInterval(t,size,verbose) ? 0 : -1); };
        bool                                       isDirty() const { return Trace<Tree>::isDirty(); };
        void                                       isDirty(bool d) const { return Trace<Tree>::isDirty(d); };
    };

}


#endif
