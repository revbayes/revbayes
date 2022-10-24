#ifndef TraceTree_H
#define TraceTree_H

#include <boost/optional.hpp>
#include <memory>
#include "Clade.h"
#include "Trace.h"
#include "Tree.h"

namespace RevBayesCore {

    class TreeSummary;

    class TraceTree : public Trace<Tree>
    {
        // These properties describe the actual trees in the trace.
        bool                                       clock = true;
        bool                                       rooted = true;

        void                                       doAddObject(Tree&& d) override;

        std::unique_ptr<TreeSummary>               the_summary;

    public:

        /*
         * Declaration of the TreeTrace class
         */
        explicit                                   TraceTree( bool c = true );
                                                   TraceTree(const TraceTree& t );
        virtual                                    ~TraceTree() = default;

        TraceTree&                                 operator=(const TraceTree&);

        const TreeSummary&                         summary() const;
        TreeSummary&                               summary();

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
