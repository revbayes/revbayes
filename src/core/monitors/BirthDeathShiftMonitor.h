#ifndef RevBayes_development_branch_BirthDeathShiftMonitor_h
#define RevBayes_development_branch_BirthDeathShiftMonitor_h

#include <iosfwd>

#include "VariableMonitor.h"


namespace RevBayesCore {
class DagNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * Monitor to print out the time spent in each state during an SSE stochastically mapped character history.
     */
    class BirthDeathShiftMonitor : public VariableMonitor {
        
    public:
        
        // Constructors and Destructors
        BirthDeathShiftMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, const std::string &del);
        BirthDeathShiftMonitor(const BirthDeathShiftMonitor &m);
        virtual ~BirthDeathShiftMonitor(void);
        
        BirthDeathShiftMonitor*                        clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                                monitorVariables(unsigned long gen);                                 //!< Monitor at generation gen
        void                                                printFileHeader(void);                                              //!< Print header
        
        // getters and setters
        void                                                swapNode(DagNode *oldN, DagNode *newN);
        
    private:
        
        // members
//        TypedDagNode<Tree>*                                 tree;
        StochasticNode<Tree>*                               bdsp;                                                              //!< The character dependent birth death process we are monitoring
    };
    
}

#endif
