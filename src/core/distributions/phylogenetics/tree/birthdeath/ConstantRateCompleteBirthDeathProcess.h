#ifndef ConstantRateCompleteBirthDeathProcess_H
#define ConstantRateCompleteBirthDeathProcess_H

#include "ConstantRateBirthDeathProcess.h"

#include <vector>

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * Constant rate complete birth-death process.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class ConstantRateCompleteBirthDeathProcess : public AbstractBirthDeathProcess {
        
    public:
        ConstantRateCompleteBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s,
                                              const TypedDagNode<double> *e, const std::string &cdt,
                                              const std::vector<Taxon> &tn, bool uo, bool mr);
        
        // public member functions
        ConstantRateCompleteBirthDeathProcess*              clone(void) const;                                                          //!< Create an independent clone
        
    protected:
        // pure virtual helper functions
        virtual double                                      computeLnProbabilityTimes(void) const;                                      //!< Compute the log-transformed probability of the current value.
        virtual double                                      lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const;         //!< Compute the log-transformed probability of the number of taxa.
        virtual double                                      lnProbTreeShape(void) const;
        virtual double                                      pSurvival(double start, double end) const;                                  //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

        // helper functions
        virtual double                                      simulateDivergenceTime(double origin, double present) const;                //!< Simulate a speciation event.
        size_t                                              num_extant( const TopologyNode & node, double * tl ) const;

        // members
        const TypedDagNode<double>*                         speciation;
        const TypedDagNode<double>*                         extinction;
        bool                                                use_mrca;
    };
    
}

#endif
