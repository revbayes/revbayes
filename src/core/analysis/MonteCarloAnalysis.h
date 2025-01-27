#ifndef MonteCarloAnalysis_H
#define MonteCarloAnalysis_H

#include "Cloneable.h"
#include "MonteCarloAnalysisOptions.h"
#include "RbFileManager.h"
#include "Parallelizable.h"
#include "RbVector.h"
#include "StoppingRule.h"
#include "Trace.h"

#include <vector>


#ifdef RB_MPI
#include <mpi.h>
#endif

namespace RevBayesCore {
    
    class Model;
    class MonteCarloSampler;
    
    /**
     * @brief Monte Carlo analysis running and managing the MonteCarloSampler objects.
     *
     * The Monte Carlo Analysis object is mostly used to run independent MonteCarloSamplers
     * and check for convergence between them.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-11
     *
     */
    class MonteCarloAnalysis : public Cloneable, public Parallelizable {
        
    public:
                
        MonteCarloAnalysis(MonteCarloSampler *m, size_t r, MonteCarloAnalysisOptions::TraceCombinationTypes ct);
        MonteCarloAnalysis(const MonteCarloAnalysis &m);
        virtual                                            ~MonteCarloAnalysis(void);                                       //!< Virtual destructor
        
        MonteCarloAnalysis&                                 operator=(const MonteCarloAnalysis &a);
        
        // public methods
        void                                                addFileMonitorExtension(const std::string &s, bool dir);
        void                                                addMonitor(const Monitor &m);
        MonteCarloAnalysis*                                 clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
#ifdef RB_MPI
        void                                                burnin(size_t g, const MPI_Comm &c, size_t ti, bool verbose=true);
#else
        void                                                burnin(size_t g, size_t ti, bool verbose=true);
#endif
        void                                                disableScreenMonitors(bool all);
        size_t                                              getCurrentGeneration(void) const;                               //!< Get the current generations number
        const Model&                                        getModel(void) const;
        void                                                initializeFromCheckpoint( const path &f );
        void                                                printPerformanceSummary(bool current_period = false) const;
        void                                                removeMonitors(void);                                           //!< Remove all monitors
#ifdef RB_MPI
        void                                                run(size_t k, RbVector<StoppingRule> r, const MPI_Comm &c, size_t ti, const path &cp_file, size_t ci=0, bool verbose=true);
#else
        void                                                run(size_t k, RbVector<StoppingRule> r, size_t ti, const path &cp_file, size_t ci=0, bool verbose=true);
#endif
#ifdef RB_MPI
        void                                                setModel(Model *m, bool redraw, const MPI_Comm &c);
#else
        void                                                setModel(Model *m, bool redraw);
#endif
        
    protected:
        void                                                setActivePIDSpecialized(size_t i, size_t n);                    //!< Set the number of processes for this class.
#ifdef RB_MPI
        void                                                resetReplicates(const MPI_Comm &c);
#else
        void                                                resetReplicates(void);
#endif

        size_t                                              replicates;
        std::vector<MonteCarloSampler*>                     runs;
        MonteCarloAnalysisOptions::TraceCombinationTypes    trace_combination;
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MonteCarloAnalysis& x);                       //!< Overloaded output operator
    
}

#endif
