#ifndef ValidationAnalysis_H
#define ValidationAnalysis_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbVector.h"

namespace RevBayesCore {
    
    class MonteCarloAnalysis;
    class Model;
    
    /**
     * @brief Posterior predictive analysis class.
     *
     * A posterior predictive analysis runs an analysis for a vector of powers
     * where the likelihood during each analysis run is raised to the given power.
     * The likelihood values and the current powers are stored in a file.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class ValidationAnalysis : public Cloneable, public Parallelizable {
        
    public:
        ValidationAnalysis(const MonteCarloAnalysis &m, size_t n);
        ValidationAnalysis(const ValidationAnalysis &a);
        virtual                                ~ValidationAnalysis(void);                               //!< Virtual destructor
        
        ValidationAnalysis&                     operator=(const ValidationAnalysis &a);
        
        // public methods
        ValidationAnalysis*                     clone(void) const;
        void                                    burnin(size_t g, size_t ti);
        void                                    runAll(size_t g);
        void                                    runSim(size_t idx, size_t g);
        void                                    summarizeAll(double c);
        void                                    summarizeSim(double c, size_t idx);
        
    private:
                
        // members
        size_t                                  num_runs;
        std::vector<MonteCarloAnalysis*>        runs;
        std::vector<Model*>                     simulation_values;

        std::map<std::string, int>              coverage_count;
    };
    
}

#endif
