#ifndef EpisodicBirthDeathProcess_H
#define EpisodicBirthDeathProcess_H

#include "BirthDeathProcess.h"

#include <vector>

namespace RevBayesCore {
    
    class Clade;
    
    class EpisodicBirthDeathProcess : public BirthDeathProcess {
        
    public:
        EpisodicBirthDeathProcess(const TypedDagNode<double> *ra,
                                  const TypedDagNode<RbVector<double> > *s,
                                  const TypedDagNode<RbVector<double> > *st,
                                  const TypedDagNode<RbVector<double> > *e,
                                  const TypedDagNode<RbVector<double> > *et,
                                  const TypedDagNode<double> *r,
                                  const TypedDagNode<double> *mp,
                                  const std::string& ss,
                                  const std::vector<Clade> &ic,
                                  const std::string &cdt,
                                  const std::vector<Taxon> &tn,
                                  Tree *t);
        
        // public member functions
        EpisodicBirthDeathProcess*                          clone(void) const  override;                                                      //!< Create an independent clone
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;                //!< Swap a parameter
        void                                                prepareProbComputation(void) const override;
        
        void                                                prepareRateIntegral(double end);                                                        //!< Compute the rate integral.
        void                                                prepareSurvivalProbability(double end, double r);                                       //!< Compute the rate integral.

        // helper functions
        double                                              lnSpeciationRate(double t) const override;
        double                                              rateIntegral(double t_low, double t_high) const override;
        double                                              computeProbabilitySurvival(double start, double end) const override;
        double                                              simulateDivergenceTime(double origin, double present) const override;                            //!< Simulate a speciation event.

    private:
        
        size_t                                              lower_index(double t) const;                                                    //!< Find the max index so that rateChangeTimes[index] < t < rateChangeTimes[index+1]
        size_t                                              lower_index(double t, size_t min, size_t max) const;                            //!< Find the max index so that rateChangeTimes[index] < t < rateChangeTimes[index+1]
        
        
        
        // members
        const TypedDagNode<RbVector<double> >*              lambda_rates;                                                                   //!< The speciation rates.
        const TypedDagNode<RbVector<double> >*              lambda_times;                                                                   //!< The time of the speciation rate changes.
        const TypedDagNode<RbVector<double> >*              mu_rates;                                                                       //!< The extinction rates.
        const TypedDagNode<RbVector<double> >*              mu_times;                                                                       //!< The times of the extinction rate changes.

        mutable std::vector<double>                         rate_change_times;
        mutable std::vector<double>                         birth;
        mutable std::vector<double>                         death;

    };
    
}

#endif
