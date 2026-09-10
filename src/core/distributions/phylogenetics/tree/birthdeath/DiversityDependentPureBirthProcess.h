#ifndef DiversityDependentPureProcess_H
#define DiversityDependentPureProcess_H

#include "AbstractBirthDeathProcess.h"

namespace RevBayesCore {
    
    class Clade;
    class Taxon;
    
    /**
     * @brief Diversity-dependent pure-birth (Yule) process.
     *
     *
     * The diversity-dependent pure-birth process is a species tree distribution without extinction events.
     * The speciation rate depends on the current diversity by the function:
     * lambda = (1.0 - n/K) * lambda_0
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class DiversityDependentPureBirthProcess : public AbstractBirthDeathProcess {
        
    public:
        DiversityDependentPureBirthProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<std::int64_t> *k,
                                           const std::string &cdt, const std::vector<Taxon> &tn);                                //!< Constructor
        
        // public member functions
        DiversityDependentPureBirthProcess*                 clone(void) const;                                                                                  //!< Create an independent clone
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
        
    private:
        
        // helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                                         //!< Compute the log-transformed probability of the current value.
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                              simulateDivergenceTime(double origin, double present) const;                                        //!< Simulate a speciation event.
        double                                              pSurvival(double start, double end) const;                                                          //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        
        // members
        const TypedDagNode<double>*                         initialSpeciation;                                                                                  //!< The initial speciation rate (lambda_0).
        const TypedDagNode<std::int64_t>*                           capacity;                                                                                           //!< The carrying capacity (K).
        
    };
    
}

#endif
