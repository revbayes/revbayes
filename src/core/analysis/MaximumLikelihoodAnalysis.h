#ifndef MaximumLikelihoodAnalysis_H
#define MaximumLikelihoodAnalysis_H

#include "Cloneable.h"
#include "MaximumLikelihoodEstimation.h"
#include "Parallelizable.h"
#include "RbVector.h"
#include "StoppingRule.h"

#include <vector>


#ifdef RB_MPI
#include <mpi.h>
#endif

namespace RevBayesCore {
    
    /**
     * @brief Monte Carlo analysis running and managing the MaximumLikelihoodEstimation objects.
     *
     * The Monte Carlo Analysis object is mostly used to run independent MaximumLikelihoodEstimations
     * and check for convergence between them.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-11
     *
     */
    class MaximumLikelihoodAnalysis : public Cloneable, public Parallelizable {
        
    public:
        MaximumLikelihoodAnalysis(MaximumLikelihoodEstimation *ms);
        MaximumLikelihoodAnalysis(const MaximumLikelihoodAnalysis &m);
        virtual                                            ~MaximumLikelihoodAnalysis(void);                                       //!< Virtual destructor
        
        MaximumLikelihoodAnalysis&                          operator=(const MaximumLikelihoodAnalysis &a);
        
        // public methods
//        void                                                addFileMonitorExtension(const std::string &s, bool dir);
//        void                                                addMonitor(const Monitor &m);
        MaximumLikelihoodAnalysis*                          clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        void                                                disableScreenMonitors(bool all);
        void                                                finishMonitors(void);
//        size_t                                              getCurrentGeneration(void) const;                               //!< Get the current generations number
        Model&                                              getModel(void);
        const Model&                                        getModel(void) const;
        void                                                monitor(size_t i) const;
//        void                                                printPerformanceSummary(void) const;
//        void                                                removeMonitors(void);                                           //!< Remove all monitors
        void                                                run(double e, bool verbose=true);
        void                                                setModel(Model *m);
        void                                                startMonitors(void);
        
    protected:
        void                                                setActivePIDSpecialized(size_t i, size_t n);                                                      //!< Set the number of processes for this class.
        
        MaximumLikelihoodEstimation*                        estimator;
        
        
#ifdef RB_MPI
        MPI_Comm analysis_comm;
#endif
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MaximumLikelihoodAnalysis& x);                       //!< Overloaded output operator
    
}

#endif
