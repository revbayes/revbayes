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
     * The per-taxon and boundary-crosser rates of Foote (2000) and the three-timer rates of
     * Alroy (2008). The value is a two by l matrix, speciation in the first row and extinction
     * in the second, with NaN wherever an interval's rate is undefined.
     */
    class FossilRateEstimatorFunction : public TypedFunction< MatrixReal > {

    public:
        enum Method { PER_TAXON, BOUNDARY_CROSSER, THREE_TIMER };                                      //!< Which of the three estimators this node computes

        FossilRateEstimatorFunction(const TypedDagNode< RbVector<Taxon> > *t,
                                    const TypedDagNode< RbVector<double> > *b,
                                    const TypedDagNode< double > *p,
                                    const TypedDagNode< double > *mx,
                                    Method m, FossilAgeAmbiguity a, bool e );
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
        bool                                    use_extant;
        mutable bool                            warned_max_age;
        mutable bool                            warned_widths;
        mutable bool                            warned_outside;
    };

}

#endif
