#ifndef RevBayes_development_branch_StochasticBranchRateMonitor_h
#define RevBayes_development_branch_StochasticBranchRateMonitor_h

#include <iosfwd>
#include <cstdint>

#include "VariableMonitor.h"


namespace RevBayesCore {
class DagNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * Monitor to print out the time spent in each state during an SSE stochastically mapped character history.
     */
    class StochasticBranchRateMonitor : public VariableMonitor {
        
    public:
        
        // Constructors and Destructors
        StochasticBranchRateMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, const std::string &del);
        StochasticBranchRateMonitor(const StochasticBranchRateMonitor &m);
        virtual ~StochasticBranchRateMonitor(void);
        
        StochasticBranchRateMonitor*                        clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                                monitorVariables(std::uint64_t gen);                                 //!< Monitor at generation gen
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
