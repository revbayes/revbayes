#ifndef BootstrapAnalysis_H
#define BootstrapAnalysis_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbVector.h"

namespace RevBayesCore {
    
    class MaximumLikelihoodAnalysis;
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
    class BootstrapAnalysis : public Cloneable, public Parallelizable {
        
    public:
        BootstrapAnalysis(const MaximumLikelihoodAnalysis &m, size_t n);
        BootstrapAnalysis(const BootstrapAnalysis &a);
        virtual                                ~BootstrapAnalysis(void);                               //!< Virtual destructor
        
        BootstrapAnalysis&                      operator=(const BootstrapAnalysis &a);
        
        // public methods
        BootstrapAnalysis*                      clone(void) const;
        void                                    runAll(double e);
        void                                    runSim(size_t idx, double e);
        
    private:
        
        // members
        size_t                                  num_runs;
        MaximumLikelihoodAnalysis*              mle_analysis;
        
    };
    
}

#endif
