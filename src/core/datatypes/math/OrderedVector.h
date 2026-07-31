#ifndef OrderedVector_H
#define OrderedVector_H

#include "RbVector.h"

#include <cstddef>
#include <vector>

namespace RevBayesCore {

    /**
     * @brief A vector of doubles that never decreases.
     *
     * The ordering is the constraint the type carries: entry j lies between its neighbours, so a
     * move that stays inside those bounds cannot break it.
     *
     * An entry is free when a move may place it anywhere its neighbours allow. The two ends are
     * not, since nothing inside the vector bounds them, and neither is an entry held equal to its
     * neighbour, where the two are one quantity with a slot to spare.
     */
    class OrderedVector {

    public:
        OrderedVector(RbVector<double> &v, const std::vector<bool> &f, const std::vector<bool> &t)
            : entries( v ), free_entry( f ), tied_to_prev( t ) {}

        size_t          size(void) const { return entries.size(); }
        double          operator[](size_t j) const { return entries[j]; }
        void            set(size_t j, double v) { entries[j] = v; }

        bool            isFree(size_t j) const
        {
            return j > 0 && j + 1 < entries.size() && j < free_entry.size() && free_entry[j];
        }

        //!< An entry held equal to j is not a bound on it, it is the same quantity, so walk past
        //!< the whole run j belongs to. Taking the adjacent entry instead gives a bound equal to
        //!< the value being moved, and the draw can then only ever go one way.
        double          lowerBound(size_t j) const
        {
            size_t k = j;
            while ( k > 0 && k < tied_to_prev.size() && tied_to_prev[k] ) --k;
            return entries[k-1];
        }
        double          upperBound(size_t j) const
        {
            size_t k = j;
            while ( k + 1 < entries.size() && k + 1 < tied_to_prev.size() && tied_to_prev[k+1] ) ++k;
            return entries[k+1];
        }

        bool            isOrdered(void) const
        {
            for (size_t j = 1; j < entries.size(); ++j) if ( entries[j] < entries[j-1] ) return false;
            return true;
        }

    private:
        RbVector<double>&           entries;
        const std::vector<bool>&    free_entry;
        const std::vector<bool>&    tied_to_prev;
    };

}

#endif
