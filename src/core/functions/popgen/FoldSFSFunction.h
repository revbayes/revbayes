#ifndef FoldSFSFunction_H
#define FoldSFSFunction_H

#include <cstdint>

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief Fold a site frequency spectrum (SFS).
     *
     * Converts an unfolded SFS of length n+1 (allele counts 0 ... n for n sampled
     * alleles) into a folded SFS of length floor(n/2)+1 (minor-allele counts
     * 0 ... floor(n/2)).
     *
     * Folding rule (0-indexed):
     *   folded[i] = unfolded[i] + unfolded[n-i]   when i != n-i
     *   folded[i] = unfolded[i]                    when i == n-i  (even n, midpoint)
     *
     * The input must have at least 2 entries (n >= 1).
     */
    class FoldSFSFunction : public TypedFunction< RbVector<std::int64_t> > {

    public:
        FoldSFSFunction(const TypedDagNode< RbVector<std::int64_t> >* unfolded_sfs);
        virtual                                       ~FoldSFSFunction(void);                                           //!< Virtual destructor

        FoldSFSFunction*                               clone(void) const;                                               //!< Create an independent clone
        void                                           update(void);

    protected:
        void                                           swapParameterInternal(const DagNode *oldP, const DagNode *newP); //!< Implementation of swapping parameters

    private:
        const TypedDagNode< RbVector<std::int64_t> >*  unfolded;
    };

}

#endif
