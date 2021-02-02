#ifndef BirthDeathBurstProcess_H
#define BirthDeathBurstProcess_H

#include <stddef.h>
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
    
    class BirthDeathBurstProcess : public AbstractBirthDeathProcess, public MemberObject< long > {
        
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
        
    protected:
        double                                              computeLnProbabilityDivergenceTimes(void);                                      //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter
        
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, long &rv) const;

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter);                                  //!< get affected nodes
        virtual void                                        keepSpecialization(DagNode* affecter);
        virtual void                                        restoreSpecialization(DagNode *restorer);
        virtual void                                        touchSpecialization(DagNode *toucher, bool touchAll);
        
        // helper functions
        double                                              computeLnProbabilityTimes(void);                                                //!< Compute the log-transformed probability of the current value.
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
