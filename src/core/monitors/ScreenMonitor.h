#ifndef ScreenMonitor_H
#define ScreenMonitor_H

#include "Monitor.h"

#include <iostream>
#include <queue>
#include <string>
#include <vector>

namespace RevBayesCore {
        
    class ScreenMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        ScreenMonitor(DagNode *n, unsigned long g, bool pp=true, bool l=true, bool pr=true);                                                                //!< Constructor with single DAG node
        ScreenMonitor(const std::vector<DagNode *> &n, unsigned long g, bool pp=true, bool l=true, bool pr=true);                                              //!< Constructor with vector of DAG node
        
        // basic methods
        ScreenMonitor*                      clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        bool                                isScreenMonitor(void) const;                                        //!< Is this a screen monitor?
        void                                monitor(unsigned long gen);                                         //!< Monitor at generation gen
        void                                reset(size_t numCycles);
        void                                setReplicateIndex(size_t idx);
        
        // ScreenMonitor functions
        void                                printHeader(void);                                                  //!< Print header

    private:
        
        // parameters
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                printWaitingTime;
        bool                                printElapsedTime;
        std::string                         prefixSeparator;
        std::string                         suffixSeparator;
        size_t                              headerPrintingInterval;
        time_t                              startTime;
        size_t                              numCycles;                                                          //!< Total number of cycles to monitor
        size_t                              currentGen;
        size_t                              startGen;
        size_t                              replicateIndex;
        
    };
    
}

#endif

