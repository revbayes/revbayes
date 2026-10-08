#ifndef FossilRateEstimatorFunction_H
#define FossilRateEstimatorFunction_H

#include <cstddef>
#include <vector>

#include "FossilIntervalCounts.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "Taxon.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore {

    /**
     * @brief Classical per-interval speciation and extinction rate estimators.
     *
     * Three methods, all reading the same per-interval taxon counts: per-taxon rates and the
     * boundary-crosser method of Foote (2000), and the three-timer method of Alroy (2008), in
     * the form Warnock et al. (2020) compare them in.
     *
     * The value is a two by l matrix, speciation in the first row and extinction in the
     * second, since both rates come from one pass over the same counts.
     *
     * A taxon's first and last appearance are the oldest and youngest intervals it was sampled
     * in, which is how these methods read a binned record.
     *
     * Intervals follow the range processes: the present, then each rate shift time, youngest
     * first. The oldest interval is unbounded unless max_age bounds it, and an unbounded
     * interval has no duration, so these estimators, which divide a proportion by one, have
     * no rate for it.
     *
     * An interval with no information returns NaN rather than throwing: the oldest interval,
     * which is unbounded and so has no width, the youngest under the three-timer method, which
     * needs both neighbours, and any interval whose denominator is empty.
     */
    class FossilRateEstimatorFunction : public TypedFunction< MatrixReal > {

    public:
        enum Method { PER_TAXON, BOUNDARY_CROSSER, THREE_TIMER };                                      //!< Which of the three estimators this node computes

        FossilRateEstimatorFunction(const TypedDagNode< RbVector<Taxon> > *t,
                                    const TypedDagNode< RbVector<double> > *b,
                                    const TypedDagNode< double > *p,
                                    const TypedDagNode< double > *mx,
                                    Method m, FossilAgeAmbiguity a);
        virtual                                ~FossilRateEstimatorFunction(void);                     //!< Virtual destructor

        FossilRateEstimatorFunction*            clone(void) const override;                            //!< Virtual copy constructor
        void                                    update(void) override;                                 //!< Recompute the rates

    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;  //!< Swap a parameter

    private:
        double                                  rate(const FossilIntervalCounts &c, size_t j, double width, double p_s, bool speciation) const;  //!< One interval's rate, NaN where undefined
        const TypedDagNode< RbVector<Taxon> >*  taxa;
        const TypedDagNode< RbVector<double> >* timeline;
        const TypedDagNode< double >*           present;
        const TypedDagNode< double >*           max_age;                                               //!< NULL leaves the oldest interval unbounded.
        Method                                  method;
        FossilAgeAmbiguity                      ambiguous;
        mutable bool                            warned_widths;
        mutable bool                            warned_outside;
    };

}

#endif
