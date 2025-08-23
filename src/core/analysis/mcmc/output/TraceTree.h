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
        bool                                       isDirty(void) const { return Trace<Tree>::isDirty(); };
        void                                       setDirty(bool d) { Trace<Tree>::setDirty(d); };
        
        int isCoveredInInterval(const std::string &v, double size, bool verbose, std::optional<bool> stochastic)
        {
            return (TreeSummary::isCoveredInInterval(v, size, stochastic.value_or(true), verbose) ? 0 : -1);            
        }; // hacky solution to make validation analyses work
    };

}


#endif
