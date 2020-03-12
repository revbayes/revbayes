#ifndef ValidationAnalysis_H
#define ValidationAnalysis_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbVector.h"

#include <map>

namespace RevBayesCore {
    
    class MonteCarloAnalysis;
    class Model;
    
    /**
     * @brief Analysis class used to run validation tests
     *
     * The validation analysis will run multiple runs of the same analysis,
     * and print the coverage of each sampled parameter (i.e. for each run,
     * whether the initial draw of the parameter is present in the credible
     * interval sampled in that run), along with the expected value of that coverage.
     *
     */
    class ValidationAnalysis : public Cloneable, public Parallelizable {
        
    public:
        ValidationAnalysis(const MonteCarloAnalysis &m, size_t n);
        ValidationAnalysis(const ValidationAnalysis &a);
        virtual                                ~ValidationAnalysis(void);
        
        ValidationAnalysis&                     operator=(const ValidationAnalysis &a);
        
        // public methods
        ValidationAnalysis*                     clone(void) const;
        void                                    burnin(size_t g, size_t ti);  //!< Perform burnin steps for all analyses.
        void                                    runAll(size_t g);  //!< Run all analyses.
        void                                    runSim(size_t idx, size_t g);  //!< Run a specific analysis.
        void                                    summarizeAll(double c);  //!< Print summary of all analyses.
        void                                    summarizeSim(double c, size_t idx);  //!< Calculate coverage counts for a specific analysis.
        
    private:
                
        // members
        size_t                                  num_runs;  //!< number of analyses to run
        std::vector<MonteCarloAnalysis*>        runs;  //!< vector of analyses
        std::vector<Model*>                     simulation_values;  //!< vector of initial values of the models

        std::map<std::string, int>              coverage_count; //!< coverage counts, indexed by parameter names
    };
    
}

#endif
