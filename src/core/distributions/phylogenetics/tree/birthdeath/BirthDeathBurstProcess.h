#ifndef BirthDeathBurstProcess_H
#define BirthDeathBurstProcess_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "MemberObject.h"
#include "RbException.h"

namespace RevBayesCore {
    
    class Taxon;
class DagNode;
template <class valueType> class RbOrderedSet;
template <class valueType> class TypedDagNode;
    
    class BirthDeathBurstProcess : public AbstractBirthDeathProcess, public MemberObject< std::int64_t > {
        
    public:
        BirthDeathBurstProcess(const TypedDagNode<double> *ra,
                               const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                               const TypedDagNode<double> *b, const TypedDagNode<double> *bt, const TypedDagNode<double> *r,
                               const std::string &cdt, const std::vector<Taxon> &tn);
        
        // public member functions
        BirthDeathBurstProcess*                             clone(void) const;
        double                                              getBurstTime(void) const;
        bool                                                isBurstSpeciation(size_t i) const;
        void                                                setBurstSpeciation(size_t i, bool tf);
        void                                                setValue(Tree *v, bool f=false);
        void                                                redrawValue(SimulationCondition c = SimulationCondition::MCMC);         //!< Draw a new random value from the distribution

    protected:
        double                                              computeLnProbabilityDivergenceTimes(void) const;                                      //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter
        
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);                                  //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // helper functions
        double                                              computeLnProbabilityTimes(void) const;                                                //!< Compute the log-transformed probability of the current value.
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                              lnProbTreeShape(void) const;
        double                                              simulateDivergenceTime(double origin, double present) const;                    //!< Simulate a speciation event.
        double                                              pSurvival(double start, double end) const;                                      //!< Compute the probability of survival of the process (without incomplete taxon sampling).
//        double                                              pZero(double t, double c1, double c2) const;
        double                                              lnQ(double t) const;
        double                                              pZero(double t) const;
        
        // members
        const TypedDagNode<double>*                         lambda;                                                                         //!< The speciation rate.
        const TypedDagNode<double>*                         mu;                                                                             //!< The extinction rate.
        const TypedDagNode<double>*                         beta;                                                                           //!< The speciation probability of at the burst time
        const TypedDagNode<double>*                         time_burst;                                                                     //!< The of the burst event
        const TypedDagNode<double>*                         rho;                                                                            //!< The sampling probability of extant taxa.
        
        std::vector<bool>                                   lineage_bursted_at_event;                                                       //!< Did the lineage burst at this event?
    };
    
    
}

#endif /* defined(BirthDeathBurstProcess_H) */
