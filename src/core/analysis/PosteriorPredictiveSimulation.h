#ifndef PosteriorPredictiveSimulation_H
#define PosteriorPredictiveSimulation_H

#include "Cloneable.h"
#include "Model.h"
#include "Parallelizable.h"
#include "RbVector.h"
#include "RbFileManager.h"
#include "Trace.h"

namespace RevBayesCore {
    
    /**
     * @brief Posterior predictive simulation from a posterior sample trace.
     *
     * This class performs a posterior predictive simulation. First, we receive a trace of 
     * all stochastic variables in the model. Then, we go through the trace and pick every i-th
     * sample. The values of the samples are used to set the values of the stochastic variables
     * in the model. You need to make sure that the trace contains values for all the stocahstic
     * nodes, otherwise the stochastic nodes will keep their current values.
     * Next, we simply simulated new values for all clamped variables. The newly simulated
     * values are written into a file per variable. We also create a directory per iteration, that is,
     * the j-th posterior predictive simulation will be written into the directory sim_j.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna and Lyndon Coghill)
     * @since Version 1.0, 2015-09-16
     *
     */
    class PosteriorPredictiveSimulation : public Cloneable, public Parallelizable {
        
    public:
        PosteriorPredictiveSimulation(const Model &m, const std::string &dir, const RbVector<ModelTrace> &t);
        PosteriorPredictiveSimulation(const Model &m, const std::string &dir, const RbVector<ModelTrace> &t, const std::vector<RevBayesCore::AncestralStateTrace> &ast);
//        PosteriorPredictiveSimulation(const PosteriorPredictiveSimulation &m);
        virtual                                            ~PosteriorPredictiveSimulation(void);                                       //!< Virtual destructor
        
//        PosteriorPredictiveSimulation&         operator=(const PosteriorPredictiveSimulation &m);                         //!< Overloaded assignment operator
        
        // public methods
        PosteriorPredictiveSimulation*                      clone(void) const;
        void                                                run(int thinning);
        
    private:
        
        Model                                               model;
        path                                                directory;
        RbVector<ModelTrace>                                traces;
        std::vector<RevBayesCore::AncestralStateTrace>      ancestral_state_traces;
        bool                                                condition_on_tips;
        
        
    };
    
}


#endif
