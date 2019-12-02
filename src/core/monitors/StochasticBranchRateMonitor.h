#ifndef RevBayes_development_branch_StochasticBranchRateMonitor_h
#define RevBayes_development_branch_StochasticBranchRateMonitor_h

#include "StateDependentSpeciationExtinctionProcess.h"
#include "VariableMonitor.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"

#include <string>
#include <vector>


namespace RevBayesCore {
    
    /**
     * Monitor to print out the time spent in each state during an SSE stochastically mapped character history.
     */
    class StochasticBranchRateMonitor : public VariableMonitor {
        
    public:
        
        // Constructors and Destructors
        StochasticBranchRateMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, const std::string &del);
        StochasticBranchRateMonitor(const StochasticBranchRateMonitor &m);
        virtual ~StochasticBranchRateMonitor(void);
        
        StochasticBranchRateMonitor*                        clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                                monitorVariables(unsigned long gen);                                 //!< Monitor at generation gen
        void                                                printFileHeader(void);                                              //!< Print header
        
        // getters and setters
        void                                                swapNode(DagNode *oldN, DagNode *newN);
        
    private:
        
        // members
//        TypedDagNode<Tree>*                                 tree;
        StochasticNode<Tree>*                               cdbdp;                                                              //!< The character dependent birth death process we are monitoring
    };
    
}

#endif
