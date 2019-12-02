#ifndef ConstantRateSerialSampledBirthDeathProcess_H
#define ConstantRateSerialSampledBirthDeathProcess_H

#include <stddef.h>
#include <iosfwd>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "RbException.h"

namespace RevBayesCore {

    class Taxon;
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;

    class ConstantRateSerialSampledBirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        ConstantRateSerialSampledBirthDeathProcess(const TypedDagNode<double> *o,
                                                const TypedDagNode<double> *s, const TypedDagNode<double> *e,
                                                const TypedDagNode<double> *p, const TypedDagNode<double> *r,
                                                const std::string &cdt, const std::vector<Taxon> &tn, bool uo,
                                                TypedDagNode<Tree> *t);

        // public member functions
        ConstantRateSerialSampledBirthDeathProcess*         clone(void) const;

    protected:
        double                                              computeLnProbabilityDivergenceTimes(void) const;                                //!< Compute the log-transformed probability of the current value.
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        // helper functions
        double                                              computeLnProbabilityTimes(void) const;                                          //!< Compute the log-transformed probability of the current value.
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                              lnProbTreeShape(void) const;
        double                                              simulateDivergenceTime(double origin, double present) const;                    //!< Simulate a speciation event.
        double                                              pSurvival(double start, double end) const;                                      //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                              pZero(double t, double c1, double c2) const;
        double                                              lnQ(double t, double c1, double c2) const;
        double                                              pHatZero(double t) const;

        // members
        const TypedDagNode<double>*                         lambda;                                                                         //!< The speciation rate.
        const TypedDagNode<double>*                         mu;                                                                             //!< The extinction rate.
        const TypedDagNode<double>*                         psi;                                                                            //!< The sampling probability of a just extinct species.
        const TypedDagNode<double>*                         rho;                                                                            //!< The sampling probability of extant taxa.

    };


}

#endif /* defined(ConstantRateSerialSampledBirthDeathProcess_H) */
