#ifndef FossilRangeMatrix_H
#define FossilRangeMatrix_H

#include "MatrixReal.h"
#include "Taxon.h"

#include <vector>

namespace RevBayesCore {

    /**
     * @brief A stratigraphic range per taxon, plus the occurrence record they were read from.
     *
     * One row per taxon, four columns ordered by age the way the timeline is, so a row reads
     * present to past and the model constraint is that it never decreases:
     *
     *     0  range_end    the range ends here (an extinction, or the marginalization limit)
     *     1  last         tau_K, the youngest augmented occurrence
     *     2  first        tau_1, the oldest augmented occurrence
     *     3  range_start  the birth, where the range separates from its ancestor
     *
     * The augmented ages are sampled, so they belong in the value rather than beside it: latent
     * state held on the distribution is dropped by a checkpoint and invisible to the DAG. The
     * occurrence record rides along as the data the ranges are conditioned on. It is clamped
     * rather than sampled, so it does not serialize -- a restored run re-reads and re-clamps it.
     */
    class FossilRangeMatrix : public MatrixReal {

    public:
        FossilRangeMatrix(void) : MatrixReal() {}
        FossilRangeMatrix(size_t n) : MatrixReal(n, 4), taxa(n) {}
        FossilRangeMatrix(const std::vector<Taxon> &t) : MatrixReal(t.size(), 4), taxa(t) {}

        FossilRangeMatrix*                  clone(void) const override { return new FossilRangeMatrix(*this); }

        double                              getRangeEnd(size_t i) const   { return (*this)[i][0]; }
        double                              getLast(size_t i) const       { return (*this)[i][1]; }
        double                              getFirst(size_t i) const      { return (*this)[i][2]; }
        double                              getRangeStart(size_t i) const { return (*this)[i][3]; }

        void                                setRangeEnd(size_t i, double v)   { (*this)[i][0] = v; }
        void                                setLast(size_t i, double v)       { (*this)[i][1] = v; }
        void                                setFirst(size_t i, double v)      { (*this)[i][2] = v; }
        void                                setRangeStart(size_t i, double v) { (*this)[i][3] = v; }

        //!< Ages never decrease along a row. A move that breaks this is repaired, not rejected.
        bool                                isOrdered(size_t i) const
        {
            return (*this)[i][0] <= (*this)[i][1] && (*this)[i][1] <= (*this)[i][2] && (*this)[i][2] <= (*this)[i][3];
        }

        const std::vector<Taxon>&           getTaxa(void) const { return taxa; }
        void                                setTaxa(const std::vector<Taxon> &t) { taxa = t; }

    private:
        std::vector<Taxon>                  taxa;
    };

}

#endif
