#ifndef ConstantRateBirthDeathProcess_H
#define ConstantRateBirthDeathProcess_H

#include "BirthDeathProcess.h"

#include <vector>

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * Constant rate birth-death process.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class ConstantRateBirthDeathProcess : public BirthDeathProcess {
        
    public:
        ConstantRateBirthDeathProcess(const TypedDagNode<double> *ra, const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                                      const TypedDagNode<double> *r, const std::string& ss, const std::vector<Clade> &ic, const std::string &cdt,
                                      const std::vector<Taxon> &tn);
        
        // public member functions
        ConstantRateBirthDeathProcess*                      clone(void) const;                                                                      //!< Create an independent clone
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        // helper functions
        double                                              lnSpeciationRate(double t) const;
        double                                              rateIntegral(double t_low, double t_high) const;
        double                                              computeProbabilitySurvival(double start, double end) const;
        virtual double                                      simulateDivergenceTime(double origin, double present) const;                //!< Simulate a speciation event.

        // members
        const TypedDagNode<double>*                         speciation;
        const TypedDagNode<double>*                         extinction;
        
    };
    
}

#endif
