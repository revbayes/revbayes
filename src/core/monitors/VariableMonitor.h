#ifndef FileMonitor_H
#define FileMonitor_H

#include <cstddef>
#include <vector>
#include <iosfwd>
#include "variant.h"

#include "AbstractFileMonitor.h"
#include "MonteCarloAnalysisOptions.h"
#include "FileFormat.h"

namespace RevBayesCore {
class DagNode;

    class VariableMonitor : public AbstractFileMonitor {

    public:
        // Constructors and Destructors
        VariableMonitor(DagNode *n, unsigned long g, const path &fname, const SampleFormat& f, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool wv=true);                                                                //!< Constructor with single DAG node
        VariableMonitor(const std::vector<DagNode *> &n, unsigned long g, const path &fname, const SampleFormat& f, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool wv=true);                                              //!< Constructor with vector of DAG node

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
        bool                                     posterior = true;
        bool                                     prior = true;
        bool                                     likelihood = true;
	std::variant<SeparatorFormat,JSONFormat> format = SeparatorFormat("\t");
    };
    
}

#endif

