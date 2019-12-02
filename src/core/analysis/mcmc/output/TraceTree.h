#ifndef TraceTree_H
#define TraceTree_H

#include "Clade.h"
#include "SBNParameters.h"
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

        bool                                       isCoveredInInterval(const std::string &v, double size, bool verbose){ return TreeSummary::isCoveredInInterval(v,size,verbose); };
        bool                                       isCoveredInInterval(const Tree &t, double size, bool verbose){ return TreeSummary::isCoveredInInterval(t,size,verbose); };
        bool                                       isDirty() const { return Trace<Tree>::isDirty(); };
        void                                       isDirty(bool d) const { return Trace<Tree>::isDirty(d); };
    };

}


#endif
