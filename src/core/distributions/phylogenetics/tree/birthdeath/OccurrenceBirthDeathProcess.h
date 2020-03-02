#ifndef OccurrenceBirthDeathProcess_H
#define OccurrenceBirthDeathProcess_H

#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "RbException.h"

namespace RevBayesCore {

class Taxon;
class DagNode;
class Tree;
class Clade;

template <class valueType> class TypedDagNode;
template <class valueType> class RbVector;

    class OccurrenceBirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        OccurrenceBirthDeathProcess(                                                      const TypedDagNode<double> *t,
                                                                                          const TypedDagNode<double> *l,
                                                                                          const TypedDagNode<double> *m,
                                                                                          const TypedDagNode<double> *p,
                                                                                          const TypedDagNode<double> *o,
                                                                                          const TypedDagNode<double> *rho,
                                                                                          const TypedDagNode<double> *r,


                                                                                          const std::string& cdt,
                                                                                          const std::vector<Taxon> &tn,
                                                                                          const TypedDagNode< RbVector<double> > *tau,
                                                                                          bool uo,

                                                                                          TypedDagNode<Tree> *tr);
virtual                                                     ~OccurrenceBirthDeathProcess(void);
        // public member functions
        OccurrenceBirthDeathProcess*                        clone(void) const;
        void                                                poolTimes(void);
        double                                              ComputeLt(void);
        double                                              ComputeMt(void);

    protected:
        double                                              computeLnProbabilityDivergenceTimes(void) ;                                //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        struct Event {
          Event(double d, std::string s, bool b) : time(d), type(s), extant(b) {};

          double time;
          std::string type;
          bool extant; // todo: get rid of this

          double getEventTime(void){ return time; }
          std::string getEventType(void){ return type; }
          bool getIsExtant(void){ return extant; }

        };

        // vector of Events
        std::vector<Event>         events;

        struct AgeCompare {
          bool operator()(const Event first, const Event second) {
                    return first.time < second.time;
          }
        };

        struct AgeCompareReverse {
          bool operator()(const Event first, const Event second) {
                    return first.time > second.time;
          }
        };


        // helper functions
        double                                              computeLnProbabilityTimes(void) const;                                          //!< Compute the log-transformed probability of the current value.
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                              lnProbTreeShape(void) const;
        double                                              simulateDivergenceTime(double origin, double present) const;                    //!< Simulate a speciation event.
        double                                              pSurvival(double start, double end) const;                                      //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        //double                                              pZero(double t, double c1, double c2) const;
        //double                                              lnQ(double t, double c1, double c2) const;
        double                                              pHatZero(double t) const;

        // members
        const TypedDagNode<double>*                         tor;                                                                            //!< Time of origin.
        const TypedDagNode<double>*                         lambda;                                                                         //!< The speciation rate.
        const TypedDagNode<double>*                         mu;                                                                             //!< The extinction rate.
        const TypedDagNode<double>*                         psi;                                                                            //!< The sampling probability of a just extinct species.
        const TypedDagNode<double>*                         omega;                                                                          //!< The probability of an observation being an occurrence.
        const TypedDagNode<double>*                         rho;                                                                            //!< The sampling probability of extant taxa.
        const TypedDagNode<double>*                         removalPr;                                                                      //!< The removal probability after sampling.
        const TypedDagNode<Tree>*                           tree;                                                                           //!< Tree
        const TypedDagNode< RbVector< double > >*           time_slices;                                                                    //!< Times at which density is computed
        int                                                 extant;                                                                         //!< Number of extant taxa
    };


}

#endif /* defined(OccurrenceBirthDeathProcess_H) */
