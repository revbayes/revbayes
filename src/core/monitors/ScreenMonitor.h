#ifndef ScreenMonitor_H
#define ScreenMonitor_H

#include <ctime>
#include <iosfwd>
#include <vector>

#include "Monitor.h"

namespace RevBayesCore {
class DagNode;

    /** @brief Monitor to write data to standard output
     *
     * @see RevLanguage::Mntr_Screen for the RL implementation
    */
    class ScreenMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        ScreenMonitor(DagNode *n, std::uint64_t g, bool pp=true, bool l=true, bool pr=true);  //!< Constructor with single DAG node
        ScreenMonitor(const std::vector<DagNode *> &n, std::uint64_t g, bool pp=true, bool l=true, bool pr=true);  //!< Constructor with vector of DAG nodes
        
        // basic methods
        ScreenMonitor*                      clone(void) const;
        
        // Monitor functions
        bool                                isScreenMonitor(void) const;
        void                                monitor(std::uint64_t gen);
        void                                reset(size_t numCycles);
        void                                setReplicateIndex(size_t idx);
        
        // ScreenMonitor functions
        void                                printHeader(void);

    private:
        
        // parameters
        bool                                posterior; //!< whether to print the posterior
        bool                                prior; //!< whether to print the prior
        bool                                likelihood; //!< whether to print the likelihood
        bool                                printWaitingTime; //!< whether to print the estimated time to completion
        bool                                printElapsedTime; //!< whether to print the time spent
        std::string                         prefixSeparator; //!< separator before each column header
        std::string                         suffixSeparator; //!<separator after each column header
        size_t                              headerPrintingInterval; //!< print the header each n iterations
        time_t                              startTime; //!< time of start of the run
        size_t                              numCycles; //!< planned number of iterations                                                       //!< Total number of cycles to monitor
        size_t                              currentGen; //!< current generation of the run
        size_t                              startGen; //!< start generation of the run
        size_t                              replicateIndex; //!< replicate index of the monitored run
        
    };
    
}

#endif

