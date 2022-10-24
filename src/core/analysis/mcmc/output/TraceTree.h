#ifndef TraceTree_H
#define TraceTree_H

#include <boost/optional.hpp>
#include "Clade.h"
#include "Trace.h"
#include "Tree.h"
#include "TreeSummary.h"

namespace RevBayesCore {

    class TraceTree : public Trace<Tree>, public TreeSummary
    {
        // These properties describe the actual trees in the trace.
        bool                                       clock = true;
        bool                                       rooted = true;

        void                                       doAddObject(Tree&& d) override;

    public:

        /*
         * Declaration of the TreeTrace class
         */
        explicit TraceTree( bool c = true );
        TraceTree(const TraceTree& t );
        virtual ~TraceTree(){}

        TraceTree*                                 clone(void) const;

        bool                                       isRooted() const;
        bool                                       isClock() const;

        int                                        isCoveredInInterval(const std::string &v, double size, bool verbose) override;
        int                                        isCoveredInInterval(const Tree &t, double size, bool verbose);

        bool                                       isDirty(void) const { return Trace<Tree>::isDirty(); };
        void                                       setDirty(bool d) { Trace<Tree>::setDirty(d); };
    };

}


#endif
