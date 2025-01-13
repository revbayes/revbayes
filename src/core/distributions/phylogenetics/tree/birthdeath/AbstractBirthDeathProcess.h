#ifndef AbstractBirthDeathProcess_H
#define AbstractBirthDeathProcess_H

#include "AbstractRootedTreeDistribution.h"

namespace RevBayesCore {
    

     /** @brief Abstract base class for birth-rate processes. */
    class AbstractBirthDeathProcess : public AbstractRootedTreeDistribution {
        
    public:
        AbstractBirthDeathProcess(const TypedDagNode<double> *ra, const std::string &cdt, const std::vector<Taxon> &tn, bool uo, Tree *t );
        
        virtual ~AbstractBirthDeathProcess(void);
        // pure virtual member functions
        virtual AbstractBirthDeathProcess*                  clone(void) const = 0;
        
    protected:
        // pure virtual helper functions
        virtual double                                      computeLnProbabilityTimes(void) const = 0;                                              //!< Compute the log-transformed probability of the current value.
        virtual size_t                                      getNumberOfTaxaAtPresent(void) const;
        virtual double                                      lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const = 0;                 //!< Compute the log-transformed probability of the number of taxa.
        virtual double                                      lnProbTreeShape(void) const;                                                            //!< Compute the log-transformed probability of the tree shape.
        virtual double                                      pSurvival(double start, double end) const = 0;                                          //!< Compute the probability of survival of the process.
        virtual void                                        prepareProbComputation(void) const;                                                     //!< Prepare likelihood calculation
        
        // helper functions
        double                                              computeLnProbabilityDivergenceTimes(void) const;                                        //!< Compute the log-transformed probability of divergence times.
        virtual std::vector<double>                         simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const;     //!< Simulate n speciation events.
        
        // members
        std::string                                         condition;                                          //!< The condition of the process (none/survival/#taxa).
    };
    
}

#endif
