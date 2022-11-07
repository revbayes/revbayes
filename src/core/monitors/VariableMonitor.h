#ifndef FileMonitor_H
#define FileMonitor_H

#include <stddef.h>
#include <vector>
#include <iosfwd>

#include "AbstractFileMonitor.h"
#include "MonteCarloAnalysisOptions.h"

namespace RevBayesCore {
class DagNode;

    class VariableMonitor : public AbstractFileMonitor {

    public:
        // Constructors and Destructors
        VariableMonitor(DagNode *n, unsigned long g, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool wv=true);                                                                //!< Constructor with single DAG node
        VariableMonitor(const std::vector<DagNode *> &n, unsigned long g, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool wv=true);                                              //!< Constructor with vector of DAG node

        // basic methods
        VariableMonitor*                        clone(void) const;                                                  //!< Clone the object
        
        // monitor methods
        virtual void                            printHeader();
        virtual void                            monitor(unsigned long gen);

        virtual void                            printFileHeader();
        virtual void                            monitorVariables(unsigned long gen);
        void                                    combineReplicates(size_t n_reps, MonteCarloAnalysisOptions::TraceCombinationTypes tc);

        // setters
        void                                    setPrintLikelihood(bool tf);
        void                                    setPrintPosterior(bool tf);
        void                                    setPrintPrior(bool tf);

    protected:
        bool                                    posterior;
        bool                                    prior;
        bool                                    likelihood;
        std::string                             separator;
    };
    
}

#endif

