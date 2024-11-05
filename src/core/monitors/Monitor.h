#ifndef Monitor_H
#define Monitor_H

#include "Cloneable.h"
#include "MonteCarloAnalysisOptions.h"

#include <vector>

namespace RevBayesCore {

    class DagNode;
    class Mcmc;
    class Model;

    /** @brief Base abstract class for all monitors.
     *
     * Monitors are tasked with saving information about one or several variable DAG node(s).
     * */

    class Monitor : public Cloneable {
    
    public:
        // Constructors and Destructors
        Monitor(std::uint64_t printgen);
        Monitor(std::uint64_t printgen, DagNode *n);
        Monitor(std::uint64_t printgen, const std::vector<DagNode *> &n);
        Monitor(const Monitor &x);
        virtual ~Monitor(void);

        // overloaded operators
        Monitor&                                    operator=(const Monitor &i);
            
        // pure virtual member functions
        virtual Monitor*                            clone(void) const = 0;
        virtual void                                monitor(std::uint64_t gen) = 0;  //!< Print one sample

        // methods you may want to overwrite
        virtual void                                addVariable(DagNode *n); //!< Add variable to monitor
        virtual void                                addFileExtension(const std::string &s, bool dir);  //!< Add extension to the file
        virtual void                                closeStream(void);  //!< Close stream after finishing writing
        virtual void                                combineReplicates(size_t n, MonteCarloAnalysisOptions::TraceCombinationTypes);  //!< Combine results from several replicate analyses
        virtual void                                disable(void);  //!< Disable this monitor
        virtual void                                enable(void);  //!< Enable this monitor
        virtual bool                                isEnabled(void) const;  //!< Is the monitor currently enabled?
        virtual bool                                isScreenMonitor(void) const;  //!< Is this a screen monitor?
        virtual bool                                isFileMonitor(void) const;  //!< Is this a file monitor?
        virtual void                                openStream(bool reopen);  //!< Open the stream for writing
        virtual void                                printHeader(void);  //!< Print header
        virtual void                                swapNode(DagNode *oldN, DagNode *newN);  //!< Replace attached node
        virtual void                                removeVariable(DagNode *n);  //!< Stop monitoring variable
        virtual void                                reset(size_t numCycles);  //!< Reset the monitored variables.

        // getters and setters
        virtual void                                setModel(Model* m);
        const std::vector<DagNode *>&               getDagNodes(void) const;
        void                                        setDagNodes(const std::vector<DagNode *>& args);
        void                                        setMcmc(Mcmc* m);

    protected:
    
        void                                        sortNodesByName(void);  //!< Sort the nodes by name
        
        // parameters
        bool                                        enabled; //!< is the monitor enabled
        std::uint64_t                               printgen;  //!< prints every printgen iterations
        Mcmc*                                       mcmc;  //!< pointer to analysis this monitor is part of
        std::vector<DagNode *>                      nodes;  //!< vector of monitored nodes
        const Model*                                model;  //!< model containing the monitored nodes
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Monitor& x);

    
}



// Global functions using the class
std::ostream&                       operator<<(std::ostream& o, const RevBayesCore::Monitor& x);   

#endif

